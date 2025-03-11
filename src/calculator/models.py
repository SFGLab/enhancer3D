from typing import Optional, List

from pydantic import BaseModel

from common.models import Enhancer3dProject, Enhancer3dProjectConfiguration


class UpsertProjectConfigurationActivityInput(BaseModel):
    project: Enhancer3dProject
    configuration: Enhancer3dProjectConfiguration


class FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(BaseModel):
    project: Enhancer3dProject
    configuration: Enhancer3dProjectConfiguration


class FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput(BaseModel):
    enhancers_promoters_chunk_paths: List[str]


class CalculateDistancesForEnhancerPromotersChunkActivityInput(BaseModel):
    project: Enhancer3dProject
    enhancers_promoters_chunk_path: str


class CalculateDistancesForEnhancerPromotersChunkActivityOutput(BaseModel):
    distances_chunk_path: str


class CalculateDistancesForProjectWorkflowInput(BaseModel):
    project: Enhancer3dProject
    configuration: Enhancer3dProjectConfiguration


class CalculateDistancesForProjectWorkflowOutput(BaseModel):
    distances_chunk_paths: List[str]
