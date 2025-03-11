import logging
import os
from concurrent.futures.thread import ThreadPoolExecutor
from typing import Tuple

import pandas as pd
from temporalio import activity

from calculator.loaders import load_enhancer_atlas_dataset_from_filesystem, load_gencode_annotation_dataset_from_filesystem
from calculator.models import FindPotentialPairsOfEnhancersPromotersForProjectActivityInput, FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput, CalculateDistancesForEnhancerPromotersChunkActivityInput, CalculateDistancesForEnhancerPromotersChunkActivityOutput
from chromatin_model.loaders.packed import load_chromatin_model_ensemble_from_filesystem
from distance_calculation.services import hydrate_enhancer_dataset_with_ensemble_data, hydrate_gencode_dataset_with_ensemble_data, extract_regional_genes_and_enhancers_for_ensemble, extract_full_genes_and_enhancers_for_ensemble, select_potential_enhances_gene_pairs, calculate_distances_for_potential_enhancer_gene_pairs
from utils.filesystem_utils import get_bucket_filesystem

logger = logging.getLogger(__name__)


@activity.defn
def find_potential_pairs_of_enhancers_promoters_for_project(input: FindPotentialPairsOfEnhancersPromotersForProjectActivityInput) -> FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput:
    bucket_fs = get_bucket_filesystem()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")
    model_repository_bucket = os.getenv("MODEL_REPOSITORY_BUCKET", "model-repository")
    logger.info(f"Using following buckets: processing={processing_bucket}, model_repository={model_repository_bucket}")

    enhancer_atlas_dataset_path = os.path.join(processing_bucket, "enhancer_atlas")
    gencode_annotation_dataset_path = os.path.join(processing_bucket, "gencode")

    project = input.project
    project_enhancer_promoter_chunks_path = os.path.join(processing_bucket, "projects", project.name, "enhancer_promoter_pairs")
    logger.info(f"Starting Enhancer3D project {project}, chunks will be output to {project_enhancer_promoter_chunks_path}")

    enhancer_atlas_dataset = load_enhancer_atlas_dataset_from_filesystem(
        bucket_fs,
        enhancer_atlas_dataset_path,
        input.enhancer_atlas_dataset_name
    )

    gencode_annotation_dataset = load_gencode_annotation_dataset_from_filesystem(
        bucket_fs,
        gencode_annotation_dataset_path,
        input.gencode_annotation_dataset_name
    )

    ensemble = load_chromatin_model_ensemble_from_filesystem(
        bucket_fs,
        model_repository_bucket,
        project.ensemble_id
    )

    # Hydrate the datasets with the ensemble data
    hydrated_enhancer_dataset = hydrate_enhancer_dataset_with_ensemble_data(
        enhancer_atlas_dataset=enhancer_atlas_dataset,
        ensemble=ensemble
    )

    hydrated_gencode_dataset = hydrate_gencode_dataset_with_ensemble_data(
        gencode_annotation_dataset=gencode_annotation_dataset,
        ensemble=ensemble
    )

    logger.info(f"Extracting genes and enhancers for the project {project}")
    regional_genes_and_enhancers_dataset = extract_regional_genes_and_enhancers_for_ensemble(
        ensemble_region=project.ensemble_region,
        ensemble=ensemble,
        hydrated_enhancer_dataset=hydrated_enhancer_dataset,
        hydrated_gencode_dataset=hydrated_gencode_dataset
    )

    logger.info(f"Extracting full genes and enhancers for the project {project}")
    full_genes_and_enhancers_dataset = extract_full_genes_and_enhancers_for_ensemble(
        hydrated_enhancer_dataset=hydrated_enhancer_dataset,
        hydrated_gencode_dataset=hydrated_gencode_dataset,
        ensemble=ensemble,
        regional_genes_and_enhancers_dataset=regional_genes_and_enhancers_dataset
    )

    logger.info(f"Joining enhancers and genes together, base_pair_linear_distance_threshold={input.base_pair_linear_distance_threshold}")
    potential_enhancer_gene_pairs = select_potential_enhances_gene_pairs(
        full_genes_and_enhancers_dataset=full_genes_and_enhancers_dataset,
        base_pair_linear_distance_threshold=input.base_pair_linear_distance_threshold
    )

    # Partition the potential enhancer-gene pairs into chunks and write them to the filesystem as binary files
    bucket_fs.makedirs(project_enhancer_promoter_chunks_path, exist_ok=True)
    chunked_potential_enhancer_gene_pairs = [
        (i, potential_enhancer_gene_pairs.iloc[i:i + input.enhancer_promoter_pairs_chunk_size])
        for i in range(0, len(potential_enhancer_gene_pairs), input.enhancer_promoter_pairs_chunk_size)
    ]

    logger.info(f"Writing enhancer-promoter pairs to the filesystem, {len(chunked_potential_enhancer_gene_pairs)} chunks")

    def write_chunk_to_filesystem(chunk_entry: Tuple[int, pd.DataFrame]) -> str:
        bucket_fs = get_bucket_filesystem()

        chunk_id, chunk = chunk_entry
        chunk_path = os.path.join(project_enhancer_promoter_chunks_path, f"chunk_{chunk_id}.parquet")

        with bucket_fs.open(chunk_path, "wb") as f:
            chunk.to_parquet(f)

        return chunk_path

    with ThreadPoolExecutor() as executor:
        enhancer_promoter_chunk_paths = list(executor.map(
            write_chunk_to_filesystem,
            chunked_potential_enhancer_gene_pairs
        ))

    return FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput(
        enhancers_promoters_chunk_paths=enhancer_promoter_chunk_paths
    )


@activity.defn
def calculate_distances_for_enhancer_promoters_chunk(input: CalculateDistancesForEnhancerPromotersChunkActivityInput) -> CalculateDistancesForEnhancerPromotersChunkActivityOutput:
    bucket_fs = get_bucket_filesystem()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")
    model_repository_bucket = os.getenv("MODEL_REPOSITORY_BUCKET", "model-repository")
    logger.info(f"Using following buckets: processing={processing_bucket}, model_repository={model_repository_bucket}")

    project = input.project
    enhancers_promoters_chunk_path = input.enhancers_promoters_chunk_path
    distances_chunk_path = os.path.join(processing_bucket, "projects", project.name, "distances", os.path.basename(enhancers_promoters_chunk_path))

    ensemble = load_chromatin_model_ensemble_from_filesystem(
        bucket_fs,
        model_repository_bucket,
        project.ensemble_id
    )

    logger.info(f"Calculating distances for enhancer-promoter pairs chunk {enhancers_promoters_chunk_path}, output will be written to {distances_chunk_path}")

    with bucket_fs.open(enhancers_promoters_chunk_path, "rb") as f:
        enhancer_promoter_pairs = pd.read_parquet(f)

    # Calculate distances for the enhancer-promoter pairs
    distances_for_potential_enhancer_gene_pairs = calculate_distances_for_potential_enhancer_gene_pairs(
        ensemble_region=project.ensemble_region,
        ensemble=ensemble,
        potential_enhancer_gene_pairs=enhancer_promoter_pairs
    )

    # Write the distances to the filesystem
    bucket_fs.makedirs(os.path.dirname(distances_chunk_path), exist_ok=True)
    with bucket_fs.open(distances_chunk_path, "wb") as f:
        distances_for_potential_enhancer_gene_pairs.to_parquet(f)

    return CalculateDistancesForEnhancerPromotersChunkActivityOutput(
        distances_chunk_path=distances_chunk_path
    )
