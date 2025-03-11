import os
from functools import lru_cache
from tempfile import TemporaryDirectory

import pandas as pd
import pyranges as pr
from fsspec import AbstractFileSystem


@lru_cache(maxsize=16)
def load_enhancer_atlas_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.csv"), "r") as dataset_file:
        return pd.read_csv(
            dataset_file,
            sep="\t",
            header=None,
            names=["chromosome", "start", "end", "score"]
        )


@lru_cache(maxsize=16)
def load_gencode_annotation_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with TemporaryDirectory() as temporary_local_directory:
        with fs.open(os.path.join(data_path, f"{dataset_name}.gtf"), "r") as dataset_file:
            # Save the dataset to a local file, necessary for pyranges to read it
            local_gencode_file_path = os.path.join(temporary_local_directory, f"{dataset_name}.gtf")
            with open(local_gencode_file_path, "w") as local_gencode_file:
                local_gencode_file.write(dataset_file.read())

            # Read the dataset using pyranges
            return pr.read_gff(
                local_gencode_file_path,
                ignore_bad=True
            )
