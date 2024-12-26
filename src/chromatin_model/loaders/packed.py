import os

import numpy as np

from chromatin_model.models import ChromatinModelEnsemble, ChromatinModelEnsembleHead


def load_chromatin_model_ensemble_from_filesystem(save_path: str, model_name: str) -> ChromatinModelEnsemble:
    model_metadata_path = os.path.join(save_path, f'{model_name}.metadata.json')
    coordinates_path = os.path.join(save_path, f'{model_name}.coordinates.npy')

    with open(model_metadata_path, 'rb') as model_metadata_fp:
        model_metadata: ChromatinModelEnsembleHead = ChromatinModelEnsembleHead.model_validate_json(model_metadata_fp.read())

    coordinates = np.load(coordinates_path)
    return ChromatinModelEnsemble.from_head(model_metadata, coordinates)
