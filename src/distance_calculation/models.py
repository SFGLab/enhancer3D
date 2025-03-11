from typing import NamedTuple

import pandas as pd


class RegionalGenesAndEnhancersDataset(NamedTuple):
    enhancers_per_ensemble_region_dataset: pd.DataFrame
    genes_per_ensemble_region_dataset: pd.DataFrame


class FullGenesAndEnhancersDataset(NamedTuple):
    full_enhancers_for_region: pd.DataFrame
    full_genes_for_region: pd.DataFrame
