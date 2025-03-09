import functools
import os
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Set

from fsspec import AbstractFileSystem
from temporalio import activity

from chromatin_model.services import repack_3dgnome_model_ensemble
from repacker.models import Repack3dgnomeModelEnsembleActivityInput, ListAllModelsInBucketActivityInput, ListAllModelsInBucketActivityOutput
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


@activity.defn
def list_all_models_in_bucket(input: ListAllModelsInBucketActivityInput) -> ListAllModelsInBucketActivityOutput:
    bucket_fs = get_bucket_filesystem()
    gnome_bucket = os.environ.get("GNOME_BUCKET", "3dgnome-landing-zone")

    # If no base paths are provided, list all models in the bucket
    base_paths_with_bucket = list(map(lambda base_path: os.path.join(gnome_bucket, base_path), input.base_paths))
    if not base_paths_with_bucket:
        base_paths_with_bucket = [gnome_bucket]

    activity.logger.info(f"Listing all models in bucket for {base_paths_with_bucket}")

    # Walk the base path to find all model paths
    def walk_model_base_path(fs: AbstractFileSystem, base_path: str) -> Set[str]:
        model_paths_under_base_path = set()
        for root, _, files in fs.walk(base_path):
            if any(file.endswith(".hcm") for file in files):
                model_paths_under_base_path.add(root)

        return model_paths_under_base_path

    # Walk the base paths in parallel to find all model paths
    model_paths = []
    with ThreadPoolExecutor() as executor:
        model_potential_paths_results = executor.map(functools.partial(walk_model_base_path, bucket_fs), base_paths_with_bucket)
        for model_potential_paths in model_potential_paths_results:
            model_paths.extend(model_potential_paths)

    # Remove the bucket prefix from the model paths
    model_paths = list(map(lambda model_path: model_path.replace(f"{gnome_bucket}/", ""), model_paths))

    activity.logger.info(f"Found {len(model_paths)} models in bucket")
    return ListAllModelsInBucketActivityOutput(model_paths=model_paths)
