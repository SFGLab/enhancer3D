from datetime import datetime
from typing import List, Sequence

from utils.pydantic_utils import BaseDatabaseModel


class DistanceCalculationEntry(BaseDatabaseModel):
    # Project
    project_id: str
    project_author: str
    project_cell_line: str
    project_executed_at: datetime

    # Ensemble
    ensemble_id: str

    # Region
    region_chr: str
    region_start: int
    region_end: int

    # Gene
    gene_id: str
    gene_chr: str
    gene_start: int
    gene_end: int
    gene_strand: str
    gene_model_position: int
    gene_model_coloring_start: int
    gene_model_coloring_end: int
    gene_TSS_pos: int

    # Enhancer
    enh_chr: str
    enh_start: int
    enh_end: int
    enh_score: float
    enh_center_position: int
    enh_model_position: int
    enh_model_coloring_start: int
    enh_model_coloring_end: int
    enh_center_pos: int
    enh_loci: str
    enh_tSS_distance: int

    # Distance
    avg_dist: float
    var_dist: float
    dist: List[float]
    number_bins: int

    @property
    def index(self) -> Sequence[str]:
        return [
            # Project identifier
            'project_id',

            # Ensemble identifier
            'ensemble_id',

            # Region identifier
            'region_chr',
            'region_start',
            'region_end',

            # Gene identifier
            'gene_id',
            'gene_start',
            'gene_end',

            # Enhancer identifier
            'enh_start',
            'enh_end'
        ]
