import os

import numpy as np

from chromatin_model import ChromatinModelEnsemble


def pack_chromatin_model_ensemble_to_filesystem(save_path: str, model: ChromatinModelEnsemble):
    if not os.path.exists(save_path):
        os.makedirs(save_path, exist_ok=True)

    model_metadata_path = os.path.join(save_path, f'{model.name}.metadata.json')
    coordinates_path = os.path.join(save_path, f'{model.name}.coordinates.npy')

    np.save(coordinates_path, model.coordinates_stack)
    with open(model_metadata_path, 'w') as model_metadata_fp:
        model_metadata_fp.write(model.head.model_dump_json(indent=4))
