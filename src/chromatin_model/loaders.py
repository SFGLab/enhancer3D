import functools
import itertools
import os
import re
from concurrent.futures import ThreadPoolExecutor
from typing import List, Iterator
from utils.pandas_utils import DataFrameBufferedSink

import numpy as np

import logging

from .models import ChromatinModel, ChromatinModelMetadata, ChromatinModelEnsemble

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


_3dgnome_model_metadata_split_pattern = re.compile(r'^(_\S+)\s+(.*)$')
_3dgnome_model_id_file_name_pattern = re.compile(r'loops_(\S+)_(\d+).hcm.smooth.cif')


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
    return int(_3dgnome_model_id_file_name_pattern.match(file_name).group(2))


def _3dgnome_model_load(model_name: str, model_id: int, model_stream: Iterator[str]) -> ChromatinModel:
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
    return ChromatinModel(
        id=model_id,
        name=model_name,
        format='cif',
        metadata=model_metadata,
        coordinates=model_coordinates
    )


def load_3dgnome_model_from_filesystem(base_path: str, model_name: str, model_id: int) -> ChromatinModel:
    if not os.path.exists(base_path):
        raise ValueError(f'Base path {base_path} does not exist')

    model_file = f'loops_{model_name}_{model_id}.hcm.smooth.cif'
    model_path = os.path.join(base_path, model_file)

    if not os.path.exists(model_path):
        raise ValueError(f'Model path {model_path} does not exist')

    with open(model_path, 'r') as model_fp:
        logger.info(f'Loading model {model_name} {model_id}')
        model_fp_stream = iter(map(str.strip, model_fp))
        return _3dgnome_model_load(
            model_name=model_name,
            model_id=model_id,
            model_stream=model_fp_stream
        )


def load_3dgnome_model_ensemble_from_filesystem(base_path: str, model_name: str) -> ChromatinModelEnsemble:
    if not os.path.exists(base_path):
        raise ValueError(f'Base path {base_path} does not exist')

    model_files = [
        f
        for f in os.listdir(base_path)
        if f.startswith(f'loops_{model_name}_') and f.endswith('.smooth.cif')
    ]

    model_ids = [
        _3dgnome_model_extract_id_from_file_name(f)
        for f in model_files
    ]

    with ThreadPoolExecutor() as executor:
        logger.info(f'Loading {len(model_ids)} models for {model_name}')
        models = list(executor.map(
            functools.partial(load_3dgnome_model_from_filesystem, base_path, model_name),
            model_ids
        ))

    metadata_stack = [m.metadata for m in models]
    coordinates_stack = (
        np
        .stack([m.coordinates for m in models])
        .astype(np.float32)
    )

    # noinspection PyTypeChecker
    return ChromatinModelEnsemble(
        name=model_name,
        format='cif',
        count=len(models),
        metadata_stack=metadata_stack,
        coordinates_stack=coordinates_stack
    )
