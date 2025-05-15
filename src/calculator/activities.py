import os
from concurrent.futures.thread import ThreadPoolExecutor
from datetime import datetime
from typing import Tuple

import pandas as pd
from temporalio import activity

from calculator.loaders import load_enhancer_atlas_dataset_from_filesystem, load_gencode_annotation_dataset_from_filesystem
from calculator.models import FindPotentialPairsOfEnhancersPromotersForProjectActivityInput, FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput, CalculateDistancesForEnhancerPromotersChunkActivityInput, CalculateDistancesForEnhancerPromotersChunkActivityOutput, UpsertProjectConfigurationActivityInput, PersistDistancesForEnhancerPromotersChunkActivityInput
from chromatin_model.loaders.packed import load_chromatin_model_ensemble_from_filesystem
from common.models import Enhancer3dProjectDatasetList
from database.models import DistanceCalculationEntry
from database.services import upsert_many_database_models_mongo, upset_many_database_models_cassandra
from distance_calculation.services import hydrate_enhancer_dataset_with_ensemble_data, hydrate_gencode_dataset_with_ensemble_data, extract_regional_genes_and_enhancers_for_ensemble, extract_full_genes_and_enhancers_for_ensemble, select_potential_enhances_gene_pairs, calculate_distances_for_potential_enhancer_gene_pairs
from utils.filesystem_utils import get_bucket_filesystem


@activity.defn
def upsert_project_configuration(input: UpsertProjectConfigurationActivityInput) -> None:
    bucket_fs = get_bucket_filesystem()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")

    project = input.project
    datasets = Enhancer3dProjectDatasetList(input.datasets)
    configuration = input.configuration

    with bucket_fs.open(os.path.join(processing_bucket, "projects", project.id, "configuration.json"), "w") as f:
        f.write(configuration.model_dump_json(indent=4))

    with bucket_fs.open(os.path.join(processing_bucket, "projects", project.id, "dataset.json"), "w") as f:
        f.write(datasets.model_dump_json(indent=4))

    with bucket_fs.open(os.path.join(processing_bucket, "projects", project.id, "project.json"), "w") as f:
        f.write(project.model_dump_json(indent=4))

    return None


@activity.defn
def find_potential_pairs_of_enhancers_promoters_for_project(input: FindPotentialPairsOfEnhancersPromotersForProjectActivityInput) -> FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput:
    bucket_fs = get_bucket_filesystem()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")
    model_repository_bucket = os.getenv("MODEL_REPOSITORY_BUCKET", "model-repository")

    enhancer_atlas_dataset_path = os.path.join(processing_bucket, "enhancer_atlas")
    gencode_annotation_dataset_path = os.path.join(processing_bucket, "gencode")

    project = input.project
    dataset = input.dataset
    configuration = input.configuration

    project_enhancer_promoter_chunks_path = os.path.join(
        processing_bucket,
        "projects",
        project.id,
        "enhancer_promoter_pairs",
        dataset.ensemble_id
    )
    activity.logger.info(f"Starting Enhancer3D project {project}, chunks will be output to {project_enhancer_promoter_chunks_path}")

    enhancer_atlas_dataset = load_enhancer_atlas_dataset_from_filesystem(
        bucket_fs,
        enhancer_atlas_dataset_path,
        dataset.enhancer_atlas_dataset_name,
        dataset.enhancer_atlas_dataset_type
    )

    gencode_annotation_dataset = load_gencode_annotation_dataset_from_filesystem(
        bucket_fs,
        gencode_annotation_dataset_path,
        dataset.gencode_annotation_dataset_name,
        dataset.gencode_annotation_dataset_type
    )

    ensemble = load_chromatin_model_ensemble_from_filesystem(
        bucket_fs,
        model_repository_bucket,
        dataset.ensemble_id
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

    activity.logger.info(f"Extracting genes and enhancers for the project {project}")
    regional_genes_and_enhancers_dataset = extract_regional_genes_and_enhancers_for_ensemble(
        ensemble_region=dataset.ensemble_region,
        ensemble=ensemble,
        hydrated_enhancer_dataset=hydrated_enhancer_dataset,
        hydrated_gencode_dataset=hydrated_gencode_dataset
    )

    activity.logger.info(f"Extracting full genes and enhancers for the project {project}")
    full_genes_and_enhancers_dataset = extract_full_genes_and_enhancers_for_ensemble(
        hydrated_enhancer_dataset=hydrated_enhancer_dataset,
        hydrated_gencode_dataset=hydrated_gencode_dataset,
        ensemble=ensemble,
        regional_genes_and_enhancers_dataset=regional_genes_and_enhancers_dataset
    )

    activity.logger.info(f"Joining enhancers and genes together, base_pair_linear_distance_threshold={configuration.base_pair_linear_distance_threshold}")
    potential_enhancer_gene_pairs = select_potential_enhances_gene_pairs(
        full_genes_and_enhancers_dataset=full_genes_and_enhancers_dataset,
        base_pair_linear_distance_threshold=configuration.base_pair_linear_distance_threshold
    )

    # Partition the potential enhancer-gene pairs into chunks and write them to the filesystem as binary files
    bucket_fs.makedirs(project_enhancer_promoter_chunks_path, exist_ok=True)
    chunked_potential_enhancer_gene_pairs = [
        (i // configuration.enhancer_promoter_pairs_chunk_size, potential_enhancer_gene_pairs[i:i + configuration.enhancer_promoter_pairs_chunk_size])
        for i in range(0, len(potential_enhancer_gene_pairs), configuration.enhancer_promoter_pairs_chunk_size)
    ]

    activity.logger.info(f"Writing enhancer-promoter pairs to the filesystem, {len(chunked_potential_enhancer_gene_pairs)} chunks")

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
        enhancers_promoters_chunk_paths=enhancer_promoter_chunk_paths,
        dataset=dataset
    )


@activity.defn
def calculate_distances_for_enhancer_promoters_chunk(input: CalculateDistancesForEnhancerPromotersChunkActivityInput) -> CalculateDistancesForEnhancerPromotersChunkActivityOutput:
    bucket_fs = get_bucket_filesystem()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")
    model_repository_bucket = os.getenv("MODEL_REPOSITORY_BUCKET", "model-repository")

    project = input.project
    dataset = input.dataset
    enhancers_promoters_chunk_path = input.enhancers_promoters_chunk_path
    distances_chunk_path = os.path.join(
        processing_bucket,
        "projects",
        project.id,
        "distances",
        dataset.ensemble_id,
        os.path.basename(enhancers_promoters_chunk_path)
    )

    ensemble = load_chromatin_model_ensemble_from_filesystem(
        bucket_fs,
        model_repository_bucket,
        dataset.ensemble_id
    )

    activity.logger.info(f"Calculating distances for enhancer-promoter pairs chunk {enhancers_promoters_chunk_path}, output will be written to {distances_chunk_path}")

    with bucket_fs.open(enhancers_promoters_chunk_path, "rb") as f:
        enhancer_promoter_pairs = pd.read_parquet(f)

    # Calculate distances for the enhancer-promoter pairs
    distances_for_potential_enhancer_gene_pairs = calculate_distances_for_potential_enhancer_gene_pairs(
        ensemble_region=dataset.ensemble_region,
        ensemble=ensemble,
        potential_enhancer_gene_pairs=enhancer_promoter_pairs
    )

    # Write the distances to the filesystem
    bucket_fs.makedirs(os.path.dirname(distances_chunk_path), exist_ok=True)
    with bucket_fs.open(distances_chunk_path, "wb") as f:
        distances_for_potential_enhancer_gene_pairs.to_parquet(f)

    return CalculateDistancesForEnhancerPromotersChunkActivityOutput(
        distances_chunk_path=distances_chunk_path,
        dataset=dataset
    )


@activity.defn
def persist_distances_for_enhancer_promoters_chunk(input: PersistDistancesForEnhancerPromotersChunkActivityInput) -> None:
    bucket_fs = get_bucket_filesystem()
    executed_at = datetime.now()

    processing_bucket = os.getenv("PROCESSING_BUCKET", "processing")
    enhancer3d_database_name = os.getenv("DATABASE_NAME", "enhancer3d")
    distance_calculation_collection_name = os.getenv("DISTANCE_CALCULATION_COLLECTION_NAME", "distance_calculation")

    project = input.project
    dataset = input.dataset
    distances_chunk_path = input.distances_chunk_path

    with bucket_fs.open(distances_chunk_path, "rb") as f:
        distances = pd.read_parquet(f)

    distances_data = [
        DistanceCalculationEntry.model_validate({
            'project_id': project.id,
            'project_authors': project.authors,
            'project_species': project.species,
            'project_cell_lines': project.cell_lines,
            'project_executed_at': executed_at,
            'ensemble_id': dataset.ensemble_id,
            **row.to_dict(),
        })
        for _, row in distances.iterrows()
    ]

    activity.logger.info(f"Persisting distances for enhancer-promoter pairs chunk {distances_chunk_path}")
    # upsert_many_database_models_mongo(
    #     database_name=enhancer3d_database_name,
    #     collection_name=distance_calculation_collection_name,
    #     data=distances_data
    # )
    upset_many_database_models_cassandra(
        keyspace=enhancer3d_database_name,
        table=distance_calculation_collection_name,
        data=distances_data
    )

    return None
