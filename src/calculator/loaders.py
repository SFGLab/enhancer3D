import os
from functools import lru_cache
from tempfile import TemporaryDirectory

import pandas as pd
import pyranges as pr
from fsspec import AbstractFileSystem

from common.models import Enhancer3dEnhancerAtlasDatasetType, Enhancer3dGencodeAnnotationDatasetType


def load_enhancer_atlas_dataset_from_filesystem_bed(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.bed"), "r") as dataset_file:
        dataset = pd.read_csv(dataset_file, sep="\t", header=None, names=["chromosome", "start", "end", "score"])
        dataset = dataset.dropna(subset=["chromosome", "start", "end", "score"])

        # enhancer id = chr:start-end
        dataset['id'] = (
            dataset['chromosome'].astype(str)
            + ':'
            + dataset['start'].astype(str)
            + '-'
            + dataset['end'].astype(str)
        )

        dataset = dataset.astype({"chromosome": "str", "start": "int", "end": "int", "score": "float"})
        return dataset


def load_enhancer_atlas_dataset_from_filesystem_tsv_liftovered_ref(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.tsv"), "r") as dataset_file:
        dataset = pd.read_csv(dataset_file, sep="\t")
        dataset = dataset.dropna(subset=["chr", "start", "end", "score"])

        dataset = dataset[["chr", "start", "end", "score"]]
        dataset = dataset.astype({"chr": "str", "start": "int", "end": "int", "score": "float"})

        # enhancer id = chr:start-end (ref)
        dataset['id'] = (
            dataset['chr'].astype(str)
            + ':'
            + dataset['start'].astype(str)
            + '-'
            + dataset['end'].astype(str)
        )

        dataset = dataset[["id", "chr", "start", "end", "score"]]
        dataset = dataset.rename(columns={"chr": "chromosome", "start": "start", "end": "end", "score": "score"})
        return dataset


def load_enhancer_atlas_dataset_from_filesystem_tsv_liftovered_mod(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.tsv"), "r") as dataset_file:
        dataset = pd.read_csv(dataset_file, sep="\t")
        dataset = dataset.dropna(subset=["chr", "start", "end", "chr_mod", "start_mod", "end_mod", "score"])

        dataset = dataset[["chr", "start", "end", "chr_mod", "start_mod", "end_mod", "score"]]
        dataset = dataset.astype({"chr": "str", "start": "int", "end": "int", "chr_mod": "str", "start_mod": "int", "end_mod": "int", "score": "float"})

        # enhancer id = chr:start-end (ref)
        dataset['id'] = (
            dataset['chr'].astype(str)
            + ':'
            + dataset['start'].astype(str)
            + '-'
            + dataset['end'].astype(str)
        )

        dataset = dataset[["id", "chr_mod", "start_mod", "end_mod", "score"]]
        dataset = dataset.rename(columns={"chr_mod": "chromosome", "start_mod": "start", "end_mod": "end", "score": "score"})
        return dataset


load_enhancer_atlas_dataset_from_filesystem_type_dispatch = {
    Enhancer3dEnhancerAtlasDatasetType.BED: load_enhancer_atlas_dataset_from_filesystem_bed,
    Enhancer3dEnhancerAtlasDatasetType.TSV_LIFTOVERED_REF: load_enhancer_atlas_dataset_from_filesystem_tsv_liftovered_ref,
    Enhancer3dEnhancerAtlasDatasetType.TSV_LIFTOVERED_MOD: load_enhancer_atlas_dataset_from_filesystem_tsv_liftovered_mod,
}


@lru_cache(maxsize=16)
def load_enhancer_atlas_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str, dataset_type: Enhancer3dEnhancerAtlasDatasetType) -> pd.DataFrame:
    return load_enhancer_atlas_dataset_from_filesystem_type_dispatch[dataset_type](fs, data_path, dataset_name)


def load_gencode_annotation_dataset_from_filesystem_gtf(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    if fs.exists(os.path.join(data_path, f"{dataset_name}-prepared.parquet")):
        with fs.open(os.path.join(data_path, f"{dataset_name}-prepared.parquet"), "rb") as prepared_dataset_file:
            return pd.read_parquet(prepared_dataset_file)

    # Load the dataset from the GTF file and prepare it for use
    with TemporaryDirectory() as temporary_local_directory:
        with fs.open(os.path.join(data_path, f"{dataset_name}.gtf"), "r") as dataset_file:
            # Save the dataset to a local file, necessary for pyranges to read it
            local_gencode_file_path = os.path.join(temporary_local_directory, f"{dataset_name}.gtf")
            with open(local_gencode_file_path, "w") as local_gencode_file:
                local_gencode_file.write(dataset_file.read())

            # Read the dataset using pyranges
            full_gencode_dataset = pr.read_gff(
                local_gencode_file_path,
                ignore_bad=True
            )

            # Prepare the dataset for use
            prepared_gencode_dataset = full_gencode_dataset[full_gencode_dataset['Feature'] == 'gene']

            # Save the prepared dataset to a parquet file
            with fs.open(os.path.join(data_path, f"{dataset_name}-prepared.parquet"), "wb") as prepared_dataset_file:
                prepared_gencode_dataset.to_parquet(prepared_dataset_file)

            return prepared_gencode_dataset[["gene_id", "Chromosome", "Start", "End", "Strand"]]


def load_gencode_annotation_dataset_from_filesystem_tsv_liftovered_ref(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.tsv"), "r") as dataset_file:
        dataset = pd.read_csv(dataset_file, sep="\t")
        dataset = dataset.dropna(subset=["gene_info", "gene_chr_ref", "gene_start_ref", "gene_end_ref", "gene_strand"])

        dataset = dataset[["gene_info", "gene_chr_ref", "gene_start_ref", "gene_end_ref", "gene_strand"]]

        # Extract gene_id from gene_info
        dataset["gene_id"] = dataset["gene_info"].str.extract(r'gene_id "([^"]+)"')[0]
        dataset = dataset.drop(columns=["gene_info"])

        dataset = dataset[["gene_id", "gene_chr_ref", "gene_start_ref", "gene_end_ref", "gene_strand"]]
        dataset = dataset.rename(columns={"gene_chr_ref": "Chromosome", "gene_start_ref": "Start", "gene_end_ref": "End", "gene_strand": "Strand"})

        # Set Feature to gene
        dataset["Feature"] = "gene"

        dataset = dataset.astype({"Chromosome": "str", "Start": "int", "End": "int", "Strand": "str"})
        return dataset


def load_gencode_annotation_dataset_from_filesystem_tsv_liftovered_mod(fs: AbstractFileSystem, data_path: str, dataset_name: str) -> pd.DataFrame:
    with fs.open(os.path.join(data_path, f"{dataset_name}.tsv"), "r") as dataset_file:
        dataset = pd.read_csv(dataset_file, sep="\t")
        dataset = dataset.dropna(subset=["gene_info", "gene_chr_mod", "gene_start_mod", "gene_end_mod", "gene_strand"])

        # Extract gene_id from gene_info
        dataset["gene_id"] = dataset["gene_info"].str.extract(r'gene_id "([^"]+)"')[0]
        dataset = dataset.drop(columns=["gene_info"])

        dataset = dataset[["gene_id", "gene_chr_mod", "gene_start_mod", "gene_end_mod", "gene_strand"]]
        dataset = dataset.rename(columns={"gene_chr_mod": "Chromosome", "gene_start_mod": "Start", "gene_end_mod": "End", "gene_strand": "Strand"})

        # Set Feature to gene
        dataset["Feature"] = "gene"

        dataset = dataset.astype({"Chromosome": "str", "Start": "int", "End": "int", "Strand": "str"})
        return dataset


load_gencode_annotation_dataset_from_filesystem_type_dispatch = {
    Enhancer3dGencodeAnnotationDatasetType.GTF: load_gencode_annotation_dataset_from_filesystem_gtf,
    Enhancer3dGencodeAnnotationDatasetType.TSV_LIFTOVERED_REF: load_gencode_annotation_dataset_from_filesystem_tsv_liftovered_ref,
    Enhancer3dGencodeAnnotationDatasetType.TSV_LIFTOVERED_MOD: load_gencode_annotation_dataset_from_filesystem_tsv_liftovered_mod,
}


@lru_cache(maxsize=16)
def load_gencode_annotation_dataset_from_filesystem(fs: AbstractFileSystem, data_path: str, dataset_name: str, dataset_type: Enhancer3dGencodeAnnotationDatasetType) -> pd.DataFrame:
    return load_gencode_annotation_dataset_from_filesystem_type_dispatch[dataset_type](fs, data_path, dataset_name)
