import pandas as pd

from chromatin_model import ChromatinModelEnsemble


def hydrate_enhancer_dataset_with_ensemble_data(
    enhancer_atlas_dataset: pd.DataFrame,
    reference_ensemble: ChromatinModelEnsemble,
    modification_ensemble: ChromatinModelEnsemble
) -> pd.DataFrame:
    # TODO: need proper algorithm for the hydration, for now assume we use the prepared files
    hydrated_enhancer_dataset = pd.read_csv(
        "./data/enhanceratlas2_liftovered_hg38_filtered_by_chrom_in_regions_GM12878_Deni_with_converted_regions.tsv",
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
    return pd.read_csv(
        "./data/gencode.v40.annotation_genes_converted_in_regions_GM12878_Deni_with_mod_regions_labelled.tsv",
        sep="\t"
    )
