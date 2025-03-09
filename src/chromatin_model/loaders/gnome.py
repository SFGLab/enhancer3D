import functools
import itertools
import os
import re
from concurrent.futures import ThreadPoolExecutor
from typing import List, Iterator, Optional, Tuple

from fsspec import AbstractFileSystem

from utils.pandas_utils import DataFrameBufferedSink

import numpy as np

import logging

from chromatin_model.models import ChromatinModel, ChromatinModelMetadata, ChromatinModelEnsemble

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


_3dgnome_model_metadata_split_pattern = re.compile(r'^(_\S+)\s+(.*)$')
_3dgnome_model_file_name_pattern = re.compile(r'loops_(\S+)_(\d+).hcm.smooth.cif')


def _3dgnome_model_extract_metadata_from_line(line: str) -> List[str]:
    # _entry.id      1
    entries = _3dgnome_model_metadata_split_pattern.match(line).groups()
    return list(map(str, entries))


def _3dgnome_model_extract_data_from_line(line: str) -> List[str]:
    # ATOM 1 C CA . ALA A 1 1 ? -7.903995 -6.946641 2.946338 1.00 99.99 C
    entries = line.split()
    return list(map(str.strip, entries))


def _3dgnome_model_extract_id_from_file_name(file_name: str) -> int:
    # loops_chr1_74280_3521135_1.hcm.smooth.cif
    return int(_3dgnome_model_file_name_pattern.match(file_name).group(2))


def _3dgnome_model_detect_name_from_path(fs: AbstractFileSystem, base_path: str) -> str:
    # find any '.smooth.cif' file
    model_file = next(
        os.path.basename(f)
        for f in fs.listdir(base_path, detail=False)
        if os.path.basename(f).endswith('.smooth.cif')
    )

    return _3dgnome_model_file_name_pattern.match(model_file).group(1)


def _3dgnome_model_load(model_name: str, model_id: int, model_stream: Iterator[str]) -> Tuple[ChromatinModelMetadata, np.ndarray]:
    model_metadata: ChromatinModelMetadata = {}
    model_columns: List[str] = []

    # Check identifier
    model_file_identifier = next(model_stream)
    if model_file_identifier != 'data_3dnome':
        raise ValueError(f'Model {model_name} is not in 3dgnome format, expected "data_3dnome" got "{model_file_identifier}"')

    # Build metadata
    for model_line in itertools.cycle(model_stream):
        if model_line == 'loop_':  # Start of columns
            break

        if model_line == '#': # Skip comments
            continue

        model_metadata_key, model_metadata_value = _3dgnome_model_extract_metadata_from_line(model_line)
        model_metadata[model_metadata_key] = model_metadata_value

    # Build columns
    first_model_line = None  # We skip this line afterward

    for model_line in model_stream:
        model_column = model_line.strip()
        if not model_column.startswith("_"):
            first_model_line = model_line
            break

        model_columns.append(model_column)

    # Build raw data
    with DataFrameBufferedSink(columns=model_columns) as sink:
        raw_model_data_entry = _3dgnome_model_extract_data_from_line(first_model_line)
        sink(raw_model_data_entry)

        for model_line in model_stream:
            raw_model_data_entry = _3dgnome_model_extract_data_from_line(model_line)
            sink(raw_model_data_entry)

        raw_model_data_df = sink.df

    model_coordinates = (
        raw_model_data_df[['_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']]
        .astype(np.float32)
        .to_numpy()
    )

    logger.info(f'Loaded model {model_name} {model_id} with {model_coordinates.shape[0]} coordinates')
    return model_metadata, model_coordinates


def _3dgnome_model_bins_load(model_name: str, model_id: int, model_bins_stream: Iterator[str]) -> Tuple[int, int, int]:
    bin_count = int(next(model_bins_stream))

    _, _, _, first_bin_value = next(model_bins_stream).split()
    first_bin_value = int(first_bin_value)

    _, _, _, second_bin_value = next(model_bins_stream).split()
    second_bin_value = int(second_bin_value)

    resolution = second_bin_value - first_bin_value

    last_bin_value = first_bin_value + (bin_count - 1) * resolution
    return first_bin_value, last_bin_value, resolution


def load_3dgnome_model_from_filesystem(fs: AbstractFileSystem, base_path: str, model_name: str, model_id: int) -> ChromatinModel:
    if not fs.exists(base_path):
        raise ValueError(f'Base path {base_path} does not exist')

    logger.info(f'Loading model {model_name} {model_id}')

    model_file = f'loops_{model_name}_{model_id}.hcm.smooth.cif'
    model_path = os.path.join(base_path, model_file)

    if not fs.exists(model_path):
        raise ValueError(f'Model path {model_path} does not exist')

    with fs.open(model_path, 'r') as model_fp:
        model_fp_stream = iter(map(str.strip, model_fp))
        model_metadata, model_coordinates = _3dgnome_model_load(
            model_name=model_name,
            model_id=model_id,
            model_stream=model_fp_stream
        )

    model_bins_file = f'loops_{model_name}_{model_id}.hcm.smooth.txt'
    model_bins_path = os.path.join(base_path, model_bins_file)

    if not fs.exists(model_bins_path):
        raise ValueError(f'Model bins path {model_bins_path} does not exist')

    with fs.open(model_bins_path, 'r') as model_bins_fp:
        model_bins_fp_stream = iter(map(str.strip, model_bins_fp))
        first_bin_value, last_bin_value, resolution = _3dgnome_model_bins_load(
            model_name=model_name,
            model_id=model_id,
            model_bins_stream=model_bins_fp_stream
        )

    return ChromatinModel(
        id=model_id,
        name=model_name,
        format='3dgnome',
        first_bin=first_bin_value,
        last_bin=last_bin_value,
        resolution=resolution,
        metadata=model_metadata,
        coordinates=model_coordinates
    )


def load_3dgnome_model_ensemble_from_filesystem(fs: AbstractFileSystem, base_path: str, model_name: Optional[str] = None) -> ChromatinModelEnsemble:
    if not fs.exists(base_path):
        raise ValueError(f'Base path {base_path} does not exist')

    if model_name is None:
        model_name = _3dgnome_model_detect_name_from_path(fs, base_path)
        logger.info(f'Detected model name {model_name} from path')

    model_ensemble_name = os.path.basename(base_path)

    model_files = [
        os.path.basename(f)
        for f in fs.listdir(base_path, detail=False)
        if os.path.basename(f).startswith(f'loops_{model_name}_') and os.path.basename(f).endswith('.smooth.cif')
    ]

    model_ids = [
        _3dgnome_model_extract_id_from_file_name(f)
        for f in model_files
    ]

    with ThreadPoolExecutor() as executor:
        logger.info(f'Loading {len(model_ids)} models for {model_name}')
        models = list(executor.map(
            functools.partial(load_3dgnome_model_from_filesystem, fs, base_path, model_name),
            model_ids
        ))

    metadata_stack = [m.metadata for m in models]
    coordinates_stack = (
        np
        .stack([m.coordinates for m in models])
        .astype(np.float32)
    )

    # noinspection PyTypeChecker
    model_ensemble = ChromatinModelEnsemble(
        name=model_name,
        format='3dgnome',
        count=len(models),
        first_bin=models[0].first_bin,
        last_bin=models[0].last_bin,
        resolution=models[0].resolution,
        metadata_stack=metadata_stack,
        coordinates_stack=coordinates_stack
    )

    logger.info(f'Loaded ensemble {model_ensemble_name} with {model_ensemble.count} models')
    return model_ensemble.rename(model_ensemble_name)
