import os
import re
from typing import List

import numpy as np
import pandas as pd

from models import ChromatinModel, ChromatinModelEnsemble, ChromatinModelMetadata

_3dgnome_model_metadata_split_pattern = re.compile(r'^(_\S+)\s+(.*)$')
_3dgnome_model_data_split_pattern = re.compile(r'\s+')


def _3dgnome_model_get_metadata_from(line: str) -> List[str]:
    return [
        str(entry)
        for entry in (
            _3dgnome_model_metadata_split_pattern
            .match(line)
            .groups()
        )
    ]


def _3dgnome_model_get_data_from(line: str) -> List[str]:
    return [
        str(entry)
        for entry in (
            _3dgnome_model_data_split_pattern
            .match(line)
            .groups()
        )
    ]


def load_3dgnome_model_from_filesystem(base_path: str, model_name: str, model_id: int) -> ChromatinModel:
    if not os.path.exists(base_path):
        raise ValueError(f'Base path {base_path} does not exist')

    model_file = f'{model_name}_{model_id}.smooth.cif'
    model_path = os.path.join(base_path, model_file)

    if not os.path.exists(model_path):
        raise ValueError(f'Model path {model_path} does not exist')

    model_metadata: ChromatinModelMetadata = {}
    model_columns: List[str] = []

    with open(model_path, 'r') as model_fp:
        # Check identifier
        model_file_identifier = model_fp.readline()
        if model_file_identifier != 'data_3dnome':
            raise ValueError(f'Model path {model_path} is not in 3dgnome format')

        # Build metadata
        for model_line in model_fp:
            if model_line == 'loop_':  # Start of columns
                break

            model_metadata_key, model_metadata_value = _3dgnome_model_get_metadata_from(model_line)
            model_metadata[model_metadata_key] = model_metadata_value

        # Build columns
        for model_line in model_fp:
            model_column = model_line.strip()
            if not model_column.startswith("_"):
                break

            model_columns.append(model_column)

        # Build raw data
        raw_model_data = pd.DataFrame(columns=model_columns)

        raw_model_data_entry = _3dgnome_model_get_data_from(model_line)
        raw_model_data = pd.concat([raw_model_data, pd.DataFrame(raw_model_data_entry)])

        for model_line in model_fp:
            raw_model_data_entry = _3dgnome_model_get_data_from(model_line)
            raw_model_data = pd.concat([raw_model_data, pd.DataFrame(raw_model_data_entry)])

    model_coordinates = (
        raw_model_data[['_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']]
        .astype(np.float32)
        .to_numpy()
    )

    return ChromatinModel(
        id=model_id,
        name=model_name,
        format='cif',
        metadata=model_metadata,
        coordinates=model_coordinates
    )
