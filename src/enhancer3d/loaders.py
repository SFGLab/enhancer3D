import os

import pandas as pd
import pyranges as pr
from fsspec import AbstractFileSystem


def load_enhancer_atlas_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.csv"), "r") as dataset_file:
        return pd.read_csv(
            dataset_file,
            sep="\t",
            header=None,
            names=["chromosome", "start", "end", "score"]
        )


def load_gencode_annotation_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.gtf"), "r") as dataset_file:
        return pr.read_gff(
            dataset_file,
            ignore_bad=True
        )
