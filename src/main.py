import argparse
import logging
import os
from typing import Tuple

import numpy as np
import pandas as pd
from pybedtools import BedTool
from typing_extensions import NamedTuple

from chromatin_model import ChromatinModelEnsemble, load_chromatin_model_ensemble_from_filesystem
from enhancer3d import ChromatinRegion, Enhancer3dProject, \
    load_enhancer_atlas_dataset_from_filesystem, load_gencode_annotation_dataset_from_filesystem, \
    hydrate_enhancer_dataset_with_ensemble_data, hydrate_gencode_dataset_with_ensemble_data

logger = logging.getLogger(__name__)

class

class RegionalGenesAndEnhancersDataset(NamedTuple):
    enhancers_per_reference_region_dataset: pd.DataFrame
    genes_per_reference_region_dataset: pd.DataFrame
    enhancers_per_modification_region_dataset: pd.DataFrame
    genes_per_modification_region_dataset: pd.DataFrame


def extract_regional_genes_and_enhancers_for_ensemble(
    reference_ensemble_region: ChromatinRegion,
    modification_ensemble_region: ChromatinRegion,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble,
    hydrated_enhancer_dataset: pd.DataFrame,
    hydrated_gencode_dataset: pd.DataFrame
) -> RegionalGenesAndEnhancersDataset:
    region_chr_ref = reference_ensemble_region.chromosome
    region_start_ref = reference_ensemble.first_bin
    region_end_ref = reference_ensemble.last_bin

    region_chr_mod = modification_ensemble_region.chromosome
    region_start_mod = modification_ensemble.first_bin
    region_end_mod = modification_ensemble.last_bin

    region_bed_ref = BedTool([(region_chr_ref, region_start_ref, region_end_ref)])
    region_bed_mod = BedTool([(region_chr_mod, region_start_mod, region_end_mod)])


    genes_bed_ref = BedTool.from_dataframe(
        hydrated_gencode_dataset
        .loc[hydrated_gencode_dataset['gene_affected_by_svs'].isna()]
        .reset_index()
        [['gene_chr_ref', 'gene_start_ref', 'gene_end_ref', 'index']]
        .astype({'gene_start_ref': 'int32', 'gene_end_ref': 'int32'})
    )

    genes_bed_mod = BedTool.from_dataframe(
        hydrated_gencode_dataset
        .loc[hydrated_gencode_dataset['gene_affected_by_svs'].isna()]
        .reset_index()
        [['gene_chr_mod', 'gene_start_mod', 'gene_end_mod', 'index']]
        .astype({'gene_start_mod': 'int32', 'gene_end_mod': 'int32'})
    )

    enhancers_bed_ref = BedTool.from_dataframe(
        hydrated_enhancer_dataset.loc[hydrated_enhancer_dataset['enh_affected_by_svs'].isna()]
        .reset_index()
        [['enh_chr_ref', 'enh_start_ref', 'enh_end_ref', 'index']]
        .astype({'enh_start_ref': 'int32', 'enh_end_ref': 'int32'})
    )

    enhancers_bed_mod = BedTool.from_dataframe(
        hydrated_enhancer_dataset.loc[hydrated_enhancer_dataset['enh_affected_by_svs'].isna()]
        .reset_index()
        [['enh_chr_mod', 'enh_start_mod', 'enh_end_mod', 'index']]
        .astype({'enh_start_mod': 'int32', 'enh_end_mod': 'int32'})
    )

    enhancers_region_ref = enhancers_bed_ref.intersect(region_bed_ref, f=1)
    genes_region_ref = genes_bed_ref.intersect(region_bed_ref, f=1)

    enhancers_region_mod = enhancers_bed_mod.intersect(region_bed_mod, f=1)
    genes_region_mod = genes_bed_mod.intersect(region_bed_mod, f=1)

    return RegionalGenesAndEnhancersDataset(
        enhancers_per_reference_region_dataset=enhancers_region_ref.to_dataframe(),
        genes_per_reference_region_dataset=genes_region_ref.to_dataframe(),
        enhancers_per_modification_region_dataset=enhancers_region_mod.to_dataframe(),
        genes_per_modification_region_dataset=genes_region_mod.to_dataframe()
    )


def extract_full_genes_and_enhancers_for_ensemble(
    hydrated_enhancer_dataset: pd.DataFrame,
    hydrated_gencode_dataset: pd.DataFrame,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble,
    regional_genes_and_enhancers_dataset: RegionalGenesAndEnhancersDataset
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    region_start_ref = reference_ensemble.first_bin
    region_start_mod = modification_ensemble.first_bin

    resolution_ref = reference_ensemble.resolution
    resolution_mod = modification_ensemble.resolution

    full_genes_for_region = (
        hydrated_gencode_dataset
        .loc[
            (hydrated_gencode_dataset.index.isin(regional_genes_and_enhancers_dataset.genes_per_reference_region_dataset['name']))
            & (hydrated_gencode_dataset.index.isin(regional_genes_and_enhancers_dataset.genes_per_modification_region_dataset['name']))
        ]
        .astype({'gene_start_mod': 'int32', 'gene_end_mod': 'int32'})
    )

    full_enhancers_for_region = (
        hydrated_enhancer_dataset
        .loc[
            (hydrated_enhancer_dataset.index.isin(regional_genes_and_enhancers_dataset.enhancers_per_reference_region_dataset['name']))
            & (hydrated_enhancer_dataset.index.isin(regional_genes_and_enhancers_dataset.enhancers_per_modification_region_dataset['name']))
        ]
        .astype({'enh_start_mod': 'int32', 'enh_end_mod': 'int32'})
    )

    full_enhancers_for_region['enh_center_position_ref'] = (full_enhancers_for_region['enh_start_ref'] + full_enhancers_for_region['enh_end_ref']) // 2
    full_enhancers_for_region['enh_center_position_mod'] = (full_enhancers_for_region['enh_start_mod'] + full_enhancers_for_region['enh_end_mod']) // 2

    full_enhancers_for_region['enh_model_position_ref'] = (full_enhancers_for_region['enh_center_position_ref'] - region_start_ref) // resolution_ref + 1
    full_enhancers_for_region['enh_model_position_mod'] = (full_enhancers_for_region['enh_center_position_mod'] - region_start_mod) // resolution_mod + 1

    full_enhancers_for_region['enh_model_coloring_start_ref'] = (full_enhancers_for_region['enh_start_ref'] - region_start_ref) // resolution_ref + 1
    full_enhancers_for_region['enh_model_coloring_end_ref'] = (full_enhancers_for_region['enh_end_ref'] - region_start_ref) // resolution_ref + 1
    full_enhancers_for_region['enh_model_coloring_start_mod'] = (full_enhancers_for_region['enh_start_mod'] - region_start_mod) // resolution_mod + 1
    full_enhancers_for_region['enh_model_coloring_end_mod'] = (full_enhancers_for_region['enh_end_mod'] - region_start_mod) // resolution_mod + 1

    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '+', 'gene_model_position_ref'] = (full_genes_for_region['gene_start_ref'] - region_start_ref) // resolution_ref + 1
    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '-', 'gene_model_position_ref'] = (full_genes_for_region['gene_end_ref'] - region_start_ref) // resolution_ref + 1

    full_genes_for_region['gene_model_coloring_start_ref'] = (full_genes_for_region['gene_start_ref'] - region_start_ref) // resolution_ref + 1
    full_genes_for_region['gene_model_coloring_end_ref'] = (full_genes_for_region['gene_end_ref'] - region_start_ref) // resolution_ref + 1

    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '+', 'gene_model_position_mod'] = (full_genes_for_region['gene_start_mod'] - region_start_mod) // resolution_mod + 1
    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '-', 'gene_model_position_mod'] = (full_genes_for_region['gene_end_mod'] - region_start_mod) // resolution_mod + 1

    full_genes_for_region['gene_model_coloring_start_mod'] = (full_genes_for_region['gene_start_mod'] - region_start_mod) // resolution_mod + 1
    full_genes_for_region['gene_model_coloring_end_mod'] = (full_genes_for_region['gene_end_mod'] - region_start_mod) // resolution_mod + 1

    full_genes_for_region = full_genes_for_region.astype({'gene_model_position_ref': 'int32', 'gene_model_position_mod': 'int32'})

    full_genes_for_region['gene_TSS_pos'] = np.where(
        full_genes_for_region['gene_strand'] == '+',
        full_genes_for_region['gene_start_ref'],
        full_genes_for_region['gene_end_ref']
    )

    full_enhancers_for_region['enh_center_pos'] = full_enhancers_for_region['enh_center_position_ref']

    full_enhancers_for_region['enh_loci'] = (
        full_enhancers_for_region['enh_chr_ref'].astype(str)
        + ":"
        + full_enhancers_for_region['enh_start_ref'].astype(str)
        + "-"
        + full_enhancers_for_region['enh_end_ref'].astype(str)
    )

    return full_enhancers_for_region, full_genes_for_region


def select_potential_enhances_gene_pairs(
    full_enhancers_for_region: pd.DataFrame,
    full_genes_for_region: pd.DataFrame
) -> pd.DataFrame:
    potential_enhancer_gene_pairs = pd.merge(
        full_genes_for_region.reset_index(),
        full_enhancers_for_region.reset_index(),
        how='cross',
        suffixes=('_gene', '_enh')
    )

    potential_enhancer_gene_pairs['gene_TSS_pos'] = np.where(
        potential_enhancer_gene_pairs['gene_strand'] == '+',
        potential_enhancer_gene_pairs['gene_start_ref'],
        potential_enhancer_gene_pairs['gene_end_ref']
    )

    potential_enhancer_gene_pairs['enh_tSS_distance'] = np.abs(
        potential_enhancer_gene_pairs['gene_TSS_pos'] - potential_enhancer_gene_pairs['enh_center_pos'])

    mask = (
        (potential_enhancer_gene_pairs['enh_tSS_distance'] <= 1000000)
        & (potential_enhancer_gene_pairs['gene_model_position_ref'] != potential_enhancer_gene_pairs['enh_model_position_ref'])
        & (potential_enhancer_gene_pairs['gene_model_position_mod'] != potential_enhancer_gene_pairs['enh_model_position_mod'])
    )

    potential_enhancer_gene_pairs = potential_enhancer_gene_pairs[mask]
    potential_enhancer_gene_pairs = potential_enhancer_gene_pairs.drop(['index_gene', 'index_enh'], axis=1)

    return potential_enhancer_gene_pairs


def calculate_distances_for_potential_enhancer_gene_pairs(
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble,
    potential_enhancer_gene_pairs: pd.DataFrame
) -> pd.DataFrame:
    # TODO: almost there
    pass


def run_distance_calculation_for_region(
    project: Enhancer3dProject,
    enhancer_atlas_dataset: pd.DataFrame,
    gencode_annotation_dataset: pd.DataFrame,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    logger.info(f"Starting Enhancer3D project: {project}")

    # Hydrate the datasets with the ensemble data
    hydrated_enhancer_dataset = hydrate_enhancer_dataset_with_ensemble_data(
        enhancer_atlas_dataset=enhancer_atlas_dataset,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble
    )

    hydrated_gencode_dataset = hydrate_gencode_dataset_with_ensemble_data(
        gencode_annotation_dataset=gencode_annotation_dataset,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble
    )

    logger.info("Extracting genes and enhancers for the reference and modification ensembles")
    regional_genes_and_enhancers_dataset = extract_regional_genes_and_enhancers_for_ensemble(
        reference_ensemble_region=project.reference_ensemble_region,
        modification_ensemble_region=project.modification_ensemble_region,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble
    )

    logger.info("Extracting full genes and enhancers for the reference and modification ensembles")
    full_enhancers_for_region, full_genes_for_region = extract_full_genes_and_enhancers_for_ensemble(
        hydrated_enhancer_dataset=hydrated_enhancer_dataset,
        hydrated_gencode_dataset=hydrated_gencode_dataset,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble,
        regional_genes_and_enhancers_dataset=regional_genes_and_enhancers_dataset
    )

    logger.info("Joining enhancers and genes together")
    potential_enhancer_gene_pairs = select_potential_enhances_gene_pairs(
        full_enhancers_for_region=full_enhancers_for_region,
        full_genes_for_region=full_genes_for_region
    )

    logger.info("Calculating distances for potential enhancer-gene pairs")
    distances_for_potential_enhancer_gene_pairs = calculate_distances_for_potential_enhancer_gene_pairs(
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble,
        potential_enhancer_gene_pairs=potential_enhancer_gene_pairs
    )

    return distances_for_potential_enhancer_gene_pairs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="enh3d")

    parser.add_argument("--project-name", type=str, required=True, default="default-project")
    parser.add_argument("--species", type=str, required=True, default="hg38")
    parser.add_argument("--cell-line", type=str, required=True, default="GM12878")

    parser.add_argument("--mapping-data-path", type=str, required=True, default="./data")
    parser.add_argument("--ensemble-data-path", type=str, required=True)

    parser.add_argument("--enhancer-atlas-dataset-name", type=str, required=True)
    parser.add_argument("--gencode-annotation-dataset-name", type=str, required=True)

    parser.add_argument("--reference-ensemble-name", type=str, required=True)
    parser.add_argument("--reference-ensemble-region", type=str, required=True)

    parser.add_argument("--modification-ensemble-name", type=str, required=True)
    parser.add_argument("--modification-ensemble-region", type=str, required=True)

    parser.add_argument("--output-path", type=str, required=True, default="./output")

    args = parser.parse_args()

    project = Enhancer3dProject(
        name=args.project_name,
        species=args.species,
        cell_line=args.cell_line,
        reference_ensemble_region=ChromatinRegion.from_string(args.reference_ensemble_region),
        modification_ensemble_region=ChromatinRegion.from_string(args.modification_ensemble_region)
    )

    mapping_data_path = args.mapping_data_path
    if not os.path.exists(mapping_data_path):
        raise FileNotFoundError(f"Mapping data path not found: {mapping_data_path}")

    enhancer_atlas_dataset = load_enhancer_atlas_dataset_from_filesystem(
        data_path=mapping_data_path,
        dataset_name=args.enhancer_atlas_dataset_name
    )

    gencode_annotation_dataset = load_gencode_annotation_dataset_from_filesystem(
        data_path=mapping_data_path,
        dataset_name=args.gencode_annotation_dataset_name
    )

    ensemble_data_path = args.ensemble_data_path
    if not os.path.exists(ensemble_data_path):
        raise FileNotFoundError(f"Reference ensemble data path not found: {ensemble_data_path}")

    reference_ensemble = load_chromatin_model_ensemble_from_filesystem(
        save_path=ensemble_data_path,
        model_name=args.reference_ensemble_name
    )

    modification_ensemble = load_chromatin_model_ensemble_from_filesystem(
        save_path=ensemble_data_path,
        model_name=args.modification_ensemble_name
    )

    reference_ensemble_region = ChromatinRegion.from_string(args.reference_ensemble_region)
    modification_ensemble_region = ChromatinRegion.from_string(args.modification_ensemble_region)

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    distances_for_potential_enhancer_gene_pairs = run_distance_calculation_for_region(
        project=project,
        enhancer_atlas_dataset=enhancer_atlas_dataset,
        gencode_annotation_dataset=gencode_annotation_dataset,
        reference_ensemble=reference_ensemble,
        modification_ensemble=modification_ensemble
    )

    logger.info(f"Saving the results to {args.output_path}")
    distances_for_potential_enhancer_gene_pairs.to_csv(
        os.path.join(args.output_path, f"{args.project_name}_distances.csv"),
        index=False
    )
