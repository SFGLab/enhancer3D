from typing import List

from pydantic import BaseModel

from common.models import Enhancer3dProject, Enhancer3dProjectConfiguration, Enhancer3dProjectDataset


class UpsertProjectConfigurationActivityInput(BaseModel):
    project: Enhancer3dProject
    datasets: List[Enhancer3dProjectDataset]
    configuration: Enhancer3dProjectConfiguration


class FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(BaseModel):
    project: Enhancer3dProject
    dataset: Enhancer3dProjectDataset
    configuration: Enhancer3dProjectConfiguration


class FindPotentialPairsOfEnhancersPromotersForProjectActivityOutput(BaseModel):
    dataset: Enhancer3dProjectDataset
    enhancers_promoters_chunk_paths: List[str]


class CalculateDistancesForEnhancerPromotersChunkActivityInput(BaseModel):
    project: Enhancer3dProject
    dataset: Enhancer3dProjectDataset
    enhancers_promoters_chunk_path: str


class CalculateDistancesForEnhancerPromotersChunkActivityOutput(BaseModel):
    dataset: Enhancer3dProjectDataset
    distances_chunk_path: str


class CalculateDistancesForProjectWorkflowInput(BaseModel):
    project: Enhancer3dProject
    datasets: List[Enhancer3dProjectDataset]
    configuration: Enhancer3dProjectConfiguration


class CalculateDistancesForProjectWorkflowOutput(BaseModel):
    distances_chunk_paths: List[str]


class PersistDistancesForEnhancerPromotersChunkActivityInput(BaseModel):
    project: Enhancer3dProject
    dataset: Enhancer3dProjectDataset
    distances_chunk_path: str


class PreloadDatasetsForProjectActivityInput(BaseModel):
    project: Enhancer3dProject
    datasets: List[Enhancer3dProjectDataset]
    configuration: Enhancer3dProjectConfiguration
