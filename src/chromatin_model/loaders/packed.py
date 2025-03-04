import os

import numpy as np
from fsspec import AbstractFileSystem

from chromatin_model.models import ChromatinModelEnsemble, ChromatinModelEnsembleHead


def load_chromatin_model_ensemble_from_filesystem(fs: AbstractFileSystem, data_path: str, model_name: str) -> ChromatinModelEnsemble:
    model_metadata_path = os.path.join(data_path, f'{model_name}.metadata.json')
    coordinates_path = os.path.join(data_path, f'{model_name}.coordinates.npy')

    with fs.open(model_metadata_path, 'rb') as model_metadata_fp:
        model_metadata: ChromatinModelEnsembleHead = ChromatinModelEnsembleHead.model_validate_json(model_metadata_fp.read())

    with fs.open(coordinates_path, 'rb') as coordinates_fp:
        coordinates = np.load(coordinates_fp)

    return ChromatinModelEnsemble.from_head(model_metadata, coordinates)
