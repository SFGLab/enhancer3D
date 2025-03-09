import logging

import pandas as pd

from chromatin_model import ChromatinModelEnsemble
from distance_calculation import Enhancer3dProject, \
    hydrate_enhancer_dataset_with_ensemble_data, hydrate_gencode_dataset_with_ensemble_data, \
    extract_regional_genes_and_enhancers_for_ensemble, extract_full_genes_and_enhancers_for_ensemble, \
    select_potential_enhances_gene_pairs, calculate_distances_for_potential_enhancer_gene_pairs


logger = logging.getLogger(__name__)


def run_distance_calculation_for_region(
    project: Enhancer3dProject,
    enhancer_atlas_dataset: pd.DataFrame,
    gencode_annotation_dataset: pd.DataFrame,
    ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    logger.info(f"Starting Enhancer3D project: {project}")

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

    logger.info("Joining enhancers and genes together")
    potential_enhancer_gene_pairs = select_potential_enhances_gene_pairs(
        full_genes_and_enhancers_dataset=full_genes_and_enhancers_dataset
    )

    logger.info("Calculating distances for potential enhancer-gene pairs")
    distances_for_potential_enhancer_gene_pairs = calculate_distances_for_potential_enhancer_gene_pairs(
        ensemble_region=project.ensemble_region,
        ensemble=ensemble,
        potential_enhancer_gene_pairs=potential_enhancer_gene_pairs
    )

    return distances_for_potential_enhancer_gene_pairs
