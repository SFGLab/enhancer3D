from typing import Optional, List

from pydantic import BaseModel


class Repack3dgnomeModelEnsembleActivityInput(BaseModel):
    source_data_path: str
    target_data_path: Optional[str] = None
    source_model_ensemble_name: Optional[str] = None
    source_model_name: Optional[str] = None


class Repack3dgnomeModelEnsembleWorkflowInput(BaseModel):
    source_data_path: str
    target_data_path: Optional[str] = None
    source_model_ensemble_name: Optional[str] = None
    source_model_name: Optional[str] = None


class ListAllModelsInBucketActivityInput(BaseModel):
    base_paths: List[str] = []


class ListAllModelsInBucketActivityOutput(BaseModel):
    model_paths: List[str]


class Repack3dgnomeManyModelEnsemblesWorkflowInput(BaseModel):
    source_paths: List[str] = []
    target_data_path: Optional[str] = None
