import logging
from typing import Optional

from fsspec import AbstractFileSystem

from chromatin_model.loaders.gnome import load_3dgnome_model_ensemble_from_filesystem
from chromatin_model.packers import pack_chromatin_model_ensemble_to_filesystem

logger = logging.getLogger(__name__)


def repack_3dgnome_model_ensemble(
    fs: AbstractFileSystem,
    source_data_path: str,
    target_data_path: str,
    source_model_name: Optional[str] = None,
) -> None:
    if not fs.exists(source_data_path):
        raise ValueError(f'Source data path {source_data_path} does not exist')

    if not fs.exists(target_data_path):
        raise ValueError(f'Target data path {target_data_path} does not exist')

    logger.info(f'Repacking 3DGNOME model ensemble from {source_data_path} to {target_data_path}')
    model_ensemble = load_3dgnome_model_ensemble_from_filesystem(fs, source_data_path, model_name=source_model_name)
    pack_chromatin_model_ensemble_to_filesystem(fs, target_data_path, model_ensemble)
