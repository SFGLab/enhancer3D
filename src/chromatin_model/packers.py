import os

import numpy as np
from fsspec import AbstractFileSystem

from chromatin_model import ChromatinModelEnsemble


def pack_chromatin_model_ensemble_to_filesystem(fs: AbstractFileSystem, save_path: str, model: ChromatinModelEnsemble):
    if not fs.exists(save_path):
        fs.makedirs(save_path, exist_ok=True)

    model_metadata_path = os.path.join(save_path, f'{model.name}.metadata.json')
    coordinates_path = os.path.join(save_path, f'{model.name}.coordinates.npy')

    with fs.open(coordinates_path, 'wb') as coordinates_fp:
        np.save(coordinates_fp, model.coordinates_stack)

    with fs.open(model_metadata_path, 'w') as model_metadata_fp:
        model_metadata_fp.write(model.head.model_dump_json(indent=4))
