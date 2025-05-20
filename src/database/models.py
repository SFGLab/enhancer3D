from datetime import datetime
from typing import List, Sequence

from common.models import Enhancer3dProject, Enhancer3dProjectDataset, Enhancer3dProjectConfiguration
from utils.pydantic_utils import BaseDatabaseModel


class ProjectConfigurationEntry(BaseDatabaseModel):
    # Project
    project_id: str

    # Configuration
    project: Enhancer3dProject
    datasets: List[Enhancer3dProjectDataset]
    configuration: Enhancer3dProjectConfiguration

    @property
    def index_fields(self) -> Sequence[str]:
        return [
            'project_id',
        ]


class DistanceCalculationEntry(BaseDatabaseModel):
    # Project
    project_id: str
    project_authors: List[str]
    project_species: List[str]
    project_cell_lines: List[str]
    project_executed_at: datetime

    # Ensemble
    ensemble_id: str

    # Region
    region_id: str
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
    gene_type: str

    # Enhancer
    enh_id: str
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
    def index_fields(self) -> Sequence[str]:
        return [
            'project_id',
            'ensemble_id',
            'region_id',
            'gene_id',
            'enh_id',
        ]
