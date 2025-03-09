from typing import Optional

from pydantic import BaseModel


class Repack3dgnomeModelEnsembleActivityInput(BaseModel):
    source_data_path: str
    target_data_path: Optional[str] = None
    source_model_name: Optional[str] = None


class Repack3dgnomeModelEnsembleWorkflowInput(BaseModel):
    source_data_path: str
    target_data_path: Optional[str] = None
    source_model_name: Optional[str] = None
