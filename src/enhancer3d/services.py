import functools
import os
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures.thread import ThreadPoolExecutor
from itertools import batched
from typing import Tuple, Dict, Any, List

import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy.stats import mannwhitneyu

from chromatin_model import ChromatinModelEnsemble
from utils.pandas_utils import DataFrameBufferedSink
from .models import ChromatinRegion, RegionalGenesAndEnhancersDataset, FullGenesAndEnhancersDataset


def hydrate_enhancer_dataset_with_ensemble_data(
    enhancer_atlas_dataset: pd.DataFrame,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    # TODO: need proper algorithm for the hydration, for now assume we use the prepared files

    data_path = os.environ.get('DATA_PATH', '.')
    hydrated_enhancer_dataset = pd.read_csv(
        os.path.join(data_path, "enhanceratlas2_liftovered_hg38_filtered_by_chrom_in_regions_GM12878_Deni_with_converted_regions.tsv"),
        sep="\t"
    )

    # Rename columns to match the expected names
    hydrated_enhancer_dataset = hydrated_enhancer_dataset.rename(columns={
        'affected_by': 'enh_affected_by_svs',
        'chr': 'enh_chr_ref',
        'start': 'enh_start_ref',
        'end': 'enh_end_ref',
        'chr_mod': 'enh_chr_mod',
        'start_mod': 'enh_start_mod',
        'end_mod': 'enh_end_mod'
    })

    return hydrated_enhancer_dataset


def hydrate_gencode_dataset_with_ensemble_data(
    gencode_annotation_dataset: pd.DataFrame,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    # TODO: need proper algorithm for the hydration, for now assume we use the prepared files

    data_path = os.environ.get('DATA_PATH', '.')
    return pd.read_csv(
        os.path.join(data_path, "gencode.v40.annotation_genes_converted_in_regions_GM12878_Deni_with_mod_regions_labelled.tsv"),
        sep="\t"
    )


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
) -> FullGenesAndEnhancersDataset:
    region_start_ref = reference_ensemble.first_bin
    region_start_mod = modification_ensemble.first_bin

    resolution_ref = reference_ensemble.resolution
    resolution_mod = modification_ensemble.resolution

    full_genes_mask = (
        (hydrated_gencode_dataset.index.isin(regional_genes_and_enhancers_dataset.genes_per_reference_region_dataset['name']))
        & (hydrated_gencode_dataset.index.isin(regional_genes_and_enhancers_dataset.genes_per_modification_region_dataset['name']))
    )
    full_genes_for_region = (
        hydrated_gencode_dataset
        .loc[full_genes_mask]
        .astype({'gene_start_mod': 'int32', 'gene_end_mod': 'int32'})
    )

    full_enhancers_mask = (
        (hydrated_enhancer_dataset.index.isin(regional_genes_and_enhancers_dataset.enhancers_per_reference_region_dataset['name']))
        & (hydrated_enhancer_dataset.index.isin(regional_genes_and_enhancers_dataset.enhancers_per_modification_region_dataset['name']))
    )
    full_enhancers_for_region = (
        hydrated_enhancer_dataset
        .loc[full_enhancers_mask]
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

    return FullGenesAndEnhancersDataset(
        full_enhancers_for_region=full_enhancers_for_region,
        full_genes_for_region=full_genes_for_region
    )


def select_potential_enhances_gene_pairs(
    full_genes_and_enhancers_dataset: FullGenesAndEnhancersDataset
) -> pd.DataFrame:
    potential_enhancer_gene_pairs = pd.merge(
        full_genes_and_enhancers_dataset.full_genes_for_region.reset_index(),
        full_genes_and_enhancers_dataset.full_enhancers_for_region.reset_index(),
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


def test_distances_for_significance(
    reference_ensemble_distances: np.ndarray,
    modification_ensemble_distances: np.ndarray,
) -> Tuple[float, float]:
    result = mannwhitneyu(reference_ensemble_distances, modification_ensemble_distances, alternative='two-sided')
    return result.pvalue, result.statistic


def calculate_distances_for_potential_enhancer_gene_pairs(
    reference_ensemble_region: ChromatinRegion,
    modification_ensemble_region: ChromatinRegion,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble,
    potential_enhancer_gene_pairs: pd.DataFrame
) -> pd.DataFrame:
    distances_dataset_sink = DataFrameBufferedSink(columns=[
        'region_chr_ref', 'region_start_ref', 'region_end_ref',
        'region_chr_mod', 'region_start_mod', 'region_end_mod',
        *potential_enhancer_gene_pairs.columns,
        'avg_dist_ref', 'avg_dist_mut', 'avg_dist_sub',
        'mwh_pvalue', 'mwh_statistic',
        'number_bins_ref', 'number_bins_mod'
    ])

    for _, potential_enhancer_gene_pair in potential_enhancer_gene_pairs.iterrows():
        reference_ensemble_distances = reference_ensemble.distance_distribution(
            potential_enhancer_gene_pair['gene_model_position_ref'],
            potential_enhancer_gene_pair['enh_model_position_ref']
        )

        modification_ensemble_distances = modification_ensemble.distance_distribution(
            potential_enhancer_gene_pair['gene_model_position_mod'],
            potential_enhancer_gene_pair['enh_model_position_mod']
        )

        mwh_pvalue, mwh_statistic = test_distances_for_significance(
            reference_ensemble_distances=reference_ensemble_distances,
            modification_ensemble_distances=modification_ensemble_distances
        )

        reference_ensemble_average_distance = reference_ensemble_distances.mean()
        modification_ensemble_average_distance = modification_ensemble_distances.mean()
        average_distance_difference = abs(reference_ensemble_average_distance - modification_ensemble_average_distance)

        distances_dataset_sink.write({
            'region_chr_ref': reference_ensemble_region.chromosome,
            'region_start_ref': reference_ensemble_region.start,
            'region_end_ref': reference_ensemble_region.end,
            'region_chr_mod': modification_ensemble_region.chromosome,
            'region_start_mod': modification_ensemble_region.start,
            'region_end_mod': modification_ensemble_region.end,
            **potential_enhancer_gene_pair.to_dict(),
            'avg_dist_ref': reference_ensemble_average_distance,
            'avg_dist_mut': modification_ensemble_average_distance,
            'avg_dist_sub': average_distance_difference,
            'mwh_pvalue': mwh_pvalue,
            'mwh_statistic': mwh_statistic,
            'number_bins_ref': reference_ensemble.count,
            'number_bins_mod': modification_ensemble.count
        })

    return distances_dataset_sink.df
