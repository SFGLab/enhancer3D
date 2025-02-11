import os

import pandas as pd
import pyranges as pr


def load_enhancer_atlas_dataset_from_filesystem(data_path: str,dataset_name: str) -> pd.DataFrame:
    return pd.read_csv(
        os.path.join(data_path, f"{dataset_name}.csv"),
        sep="\t",
        header=None,
        names=["chromosome", "start", "end", "score"]
    )


def load_gencode_annotation_dataset_from_filesystem(data_path: str,dataset_name: str) -> pd.DataFrame:
    return pr.read_gff(
        f=os.path.join(data_path, f"{dataset_name}.gtf"),
        as_df=True,
        ignore_bad=True
    )
