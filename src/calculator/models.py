from typing import Optional, List

from pydantic import BaseModel

from common.models import Enhancer3dProject


class FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(BaseModel):
    project: Enhancer3dProject
    enhancer_atlas_dataset_name: str
    gencode_annotation_dataset_name: str
    base_pair_linear_distance_threshold: Optional[int] = None
    enhancer_promoter_pairs_chunk_size: int = 1000


class FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput(BaseModel):
    enhancers_promoters_chunk_paths: List[str]


class CalculateDistancesForEnhancerPromotersChunkActivityInput(BaseModel):
    project: Enhancer3dProject
    enhancers_promoters_chunk_path: str


class CalculateDistancesForEnhancerPromotersChunkActivityOutput(BaseModel):
    distances_chunk_path: str


class CalculateDistancesForProjectWorkflowInput(BaseModel):
    project: Enhancer3dProject
    enhancer_atlas_dataset_name: str
    gencode_annotation_dataset_name: str
    base_pair_linear_distance_threshold: Optional[int] = None
    enhancer_promoter_pairs_chunk_size: int = 1000


class CalculateDistancesForProjectWorkflowOutput(BaseModel):
    distances_chunk_paths: List[str]
