import os

from temporalio import activity

from chromatin_model.services import repack_3dgnome_model_ensemble
from repacker.models import Repack3dgnomeModelEnsembleActivityInput
from utils.filesystem_utils import get_bucket_filesystem


@activity.defn
def repack_3dgnome_model_ensemble_from_bucket(input: Repack3dgnomeModelEnsembleActivityInput) -> None:
    bucket_fs = get_bucket_filesystem()
    gnome_bucket = os.environ.get("GNOME_BUCKET", "3dgnome-landing-zone")
    model_repository_bucket = os.environ.get("MODEL_REPOSITORY_BUCKET", "model-repository")

    source_data_path_with_bucket = os.path.join(gnome_bucket, input.source_data_path)
    target_data_path_with_bucket = (
        os.path.join(model_repository_bucket, input.source_model_name)
        if input.target_data_path
        else model_repository_bucket
    )

    activity.logger.info(f"Repacking 3D-GNOME model ensemble from {source_data_path_with_bucket} to {target_data_path_with_bucket}")
    repack_3dgnome_model_ensemble(
        fs=bucket_fs,
        source_data_path=source_data_path_with_bucket,
        target_data_path=target_data_path_with_bucket,
        source_model_name=input.source_model_name,
    )

    activity.logger.info(f"Repacking complete")
