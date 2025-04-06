from typing import Optional

import numpy as np
import pandas as pd
import pyranges as pr

from chromatin_model.models import ChromatinModelEnsemble
from common.models import ChromatinRegion
from utils.pandas_utils import DataFrameBufferedSink
from distance_calculation.models import RegionalGenesAndEnhancersDataset, FullGenesAndEnhancersDataset


def hydrate_enhancer_dataset_with_ensemble_data(
    enhancer_atlas_dataset: pd.DataFrame,
    ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    hydrated_enhancer_dataset = enhancer_atlas_dataset.rename(columns={
        'id': 'enh_id',
        'chromosome': 'enh_chr',
        'start': 'enh_start',
        'end': 'enh_end',
        'score': 'enh_score',
    })

    hydrated_enhancer_dataset = hydrated_enhancer_dataset[['enh_id', 'enh_chr', 'enh_start', 'enh_end', 'enh_score']]
    return hydrated_enhancer_dataset


def hydrate_gencode_dataset_with_ensemble_data(
    gencode_annotation_dataset: pd.DataFrame,
    ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    hydrated_gencode_dataset = gencode_annotation_dataset[gencode_annotation_dataset['Feature'] == 'gene']
    hydrated_gencode_dataset = hydrated_gencode_dataset.rename(columns={
        'gene_id': 'gene_id',
        'Chromosome': 'gene_chr',
        'Start': 'gene_start',
        'End': 'gene_end',
        'Strand': 'gene_strand',
    })

    hydrated_gencode_dataset = hydrated_gencode_dataset[['gene_id', 'gene_chr', 'gene_start', 'gene_end', 'gene_strand']]
    return hydrated_gencode_dataset


def extract_regional_genes_and_enhancers_for_ensemble(
    ensemble_region: ChromatinRegion,
    ensemble: ChromatinModelEnsemble,
    hydrated_enhancer_dataset: pd.DataFrame,
    hydrated_gencode_dataset: pd.DataFrame
) -> RegionalGenesAndEnhancersDataset:
    region_chr = ensemble_region.chromosome
    region_start = ensemble.first_bin
    region_end = ensemble.last_bin

    region_pr = pr.PyRanges({
        "Chromosome": [region_chr],
        "Start": [region_start],
        "End": [region_end],
    })

    enhancers_pr = pr.PyRanges(
        hydrated_enhancer_dataset.rename(columns={
            "enh_chr": "Chromosome",
            "enh_start": "Start",
            "enh_end": "End",
        })
    )

    genes_pr = pr.PyRanges(
        hydrated_gencode_dataset.rename(columns={
            "gene_chr": "Chromosome",
            "gene_start": "Start",
            "gene_end": "End",
        })
    )

    genes_region = genes_pr.intersect(region_pr, multiple="containment")
    enhancers_region = enhancers_pr.intersect(region_pr, multiple="containment")

    return RegionalGenesAndEnhancersDataset(
        enhancers_per_ensemble_region_dataset=enhancers_region.rename({
            "Chromosome": "enh_chr",
            "Start": "enh_start",
            "End": "enh_end",
        }),
        genes_per_ensemble_region_dataset=genes_region.rename({
            "Chromosome": "gene_chr",
            "Start": "gene_start",
            "End": "gene_end",
        })
    )


def extract_full_genes_and_enhancers_for_ensemble(
    hydrated_enhancer_dataset: pd.DataFrame,
    hydrated_gencode_dataset: pd.DataFrame,
    ensemble: ChromatinModelEnsemble,
    regional_genes_and_enhancers_dataset: RegionalGenesAndEnhancersDataset
) -> FullGenesAndEnhancersDataset:
    region_start = ensemble.first_bin
    resolution = ensemble.resolution

    full_genes_mask = hydrated_gencode_dataset.index.isin(regional_genes_and_enhancers_dataset.genes_per_ensemble_region_dataset.index)
    full_genes_for_region = (
        hydrated_gencode_dataset
        .loc[full_genes_mask]
        .astype({'gene_start': 'int32', 'gene_end': 'int32'})
    )

    full_enhancers_mask = hydrated_enhancer_dataset.index.isin(regional_genes_and_enhancers_dataset.enhancers_per_ensemble_region_dataset.index)
    full_enhancers_for_region = (
        hydrated_enhancer_dataset
        .loc[full_enhancers_mask]
        .astype({'enh_start': 'int32', 'enh_end': 'int32'})
    )

    full_enhancers_for_region['enh_center_position'] = (full_enhancers_for_region['enh_start'] + full_enhancers_for_region['enh_end']) // 2
    full_enhancers_for_region['enh_model_position'] = (full_enhancers_for_region['enh_center_position'] - region_start) // resolution + 1

    full_enhancers_for_region['enh_model_coloring_start'] = (full_enhancers_for_region['enh_start'] - region_start) // resolution + 1
    full_enhancers_for_region['enh_model_coloring_end'] = (full_enhancers_for_region['enh_end'] - region_start) // resolution + 1

    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '+', 'gene_model_position'] = (full_genes_for_region['gene_start'] - region_start) // resolution + 1
    full_genes_for_region.loc[full_genes_for_region['gene_strand'] == '-', 'gene_model_position'] = (full_genes_for_region['gene_end'] - region_start) // resolution + 1

    full_genes_for_region['gene_model_coloring_start'] = (full_genes_for_region['gene_start'] - region_start) // resolution + 1
    full_genes_for_region['gene_model_coloring_end'] = (full_genes_for_region['gene_end'] - region_start) // resolution + 1

    full_genes_for_region = full_genes_for_region.astype({'gene_model_position': 'int32'})

    full_genes_for_region['gene_TSS_pos'] = np.where(
        full_genes_for_region['gene_strand'] == '+',
        full_genes_for_region['gene_start'],
        full_genes_for_region['gene_end']
    )

    full_enhancers_for_region['enh_center_pos'] = full_enhancers_for_region['enh_center_position']

    full_enhancers_for_region['enh_loci'] = (
        full_enhancers_for_region['enh_chr'].astype(str)
        + ":"
        + full_enhancers_for_region['enh_start'].astype(str)
        + "-"
        + full_enhancers_for_region['enh_end'].astype(str)
    )

    return FullGenesAndEnhancersDataset(
        full_enhancers_for_region=full_enhancers_for_region,
        full_genes_for_region=full_genes_for_region
    )


def select_potential_enhances_gene_pairs(
    full_genes_and_enhancers_dataset: FullGenesAndEnhancersDataset,
    base_pair_linear_distance_threshold: Optional[int] = None
) -> pd.DataFrame:
    potential_enhancer_gene_pairs = pd.merge(
        full_genes_and_enhancers_dataset.full_genes_for_region.reset_index(),
        full_genes_and_enhancers_dataset.full_enhancers_for_region.reset_index(),
        how='cross',
        suffixes=('_gene', '_enh')
    )

    potential_enhancer_gene_pairs['gene_TSS_pos'] = np.where(
        potential_enhancer_gene_pairs['gene_strand'] == '+',
        potential_enhancer_gene_pairs['gene_start'],
        potential_enhancer_gene_pairs['gene_end']
    )

    potential_enhancer_gene_pairs['enh_tSS_distance'] = np.abs(
        potential_enhancer_gene_pairs['gene_TSS_pos'] - potential_enhancer_gene_pairs['enh_center_pos'])

    mask = (potential_enhancer_gene_pairs['gene_model_position'] != potential_enhancer_gene_pairs['enh_model_position'])
    if base_pair_linear_distance_threshold is not None:
        mask = mask & (potential_enhancer_gene_pairs['enh_tSS_distance'] <= base_pair_linear_distance_threshold)

    potential_enhancer_gene_pairs = potential_enhancer_gene_pairs[mask]
    potential_enhancer_gene_pairs = potential_enhancer_gene_pairs.drop(['index_gene', 'index_enh'], axis=1)

    return potential_enhancer_gene_pairs


def calculate_distances_for_potential_enhancer_gene_pairs(
    ensemble_region: ChromatinRegion,
    ensemble: ChromatinModelEnsemble,
    potential_enhancer_gene_pairs: pd.DataFrame
) -> pd.DataFrame:
    distances_dataset_sink = DataFrameBufferedSink(columns=[
        'region_id', 'region_chr', 'region_start', 'region_end',
        *potential_enhancer_gene_pairs.columns,
        'avg_dist', 'var_dist', 'dist', 'number_bins'
    ])

    for _, potential_enhancer_gene_pair in potential_enhancer_gene_pairs.iterrows():
        ensemble_distances = ensemble.distance_distribution(
            potential_enhancer_gene_pair['gene_model_position'],
            potential_enhancer_gene_pair['enh_model_position']
        )

        ensemble_average_distance = ensemble_distances.mean()
        ensemble_distance_variance = ensemble_distances.var()

        distances_dataset_sink.write({
            'region_id': str(ensemble_region),
            'region_chr': ensemble_region.chromosome,
            'region_start': ensemble_region.start,
            'region_end': ensemble_region.end,
            **potential_enhancer_gene_pair.to_dict(),
            'avg_dist': ensemble_average_distance,
            'var_dist': ensemble_distance_variance,
            'dist': ensemble_distances,
            'number_bins': ensemble.count
        })

    return distances_dataset_sink.df
