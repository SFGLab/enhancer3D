from typing import Dict, Any, List

import numpy as np
from pydantic import BaseModel
from numpydantic import NDArray, Shape

ChromatinModelMetadata = Dict[str, Any]


class ChromatinModel(BaseModel):
    id: int
    name: str
    format: str
    first_bin: int
    last_bin: int
    resolution: int
    metadata: ChromatinModelMetadata
    coordinates: NDArray[Shape["*, 3"], np.float32]


class ChromatinModelEnsembleHead(BaseModel):
    name: str
    format: str
    first_bin: int
    last_bin: int
    resolution: int
    count: int

    metadata_stack: List[ChromatinModelMetadata]


class ChromatinModelEnsemble(BaseModel):
    name: str
    format: str
    first_bin: int
    last_bin: int
    resolution: int
    count: int

    metadata_stack: List[ChromatinModelMetadata]
    coordinates_stack: NDArray[Shape["*, *, 3"], np.float32]

    @property
    def head(self) -> ChromatinModelEnsembleHead:
        return ChromatinModelEnsembleHead(
            name=self.name,
            format=self.format,
            first_bin=self.first_bin,
            last_bin=self.last_bin,
            resolution=self.resolution,
            count=self.count,
            metadata_stack=self.metadata_stack
        )

    def rename(self, new_name: str) -> 'ChromatinModelEnsemble':
        return self.model_copy(update={'name': new_name})

    def get_model(self, model_id: int) -> ChromatinModel:
        if not 0 <= model_id < self.count:
            raise ValueError(f'Model ensemble does not have {model_id} model.')

        coordinates = self.coordinates_stack[model_id]
        metadata = self.metadata_stack[model_id]
        return ChromatinModel(
            id=model_id,
            name=self.name,
            format=self.format,
            first_bin=self.first_bin,
            last_bin=self.last_bin,
            resolution=self.resolution,
            metadata=metadata,
            coordinates=coordinates
        )

    @staticmethod
    def from_models(models: List[ChromatinModel]) -> 'ChromatinModelEnsemble':
        metadata_stack = [m.metadata for m in models]
        coordinates_stack = np.stack([m.coordinates for m in models])

        # noinspection PyTypeChecker
        return ChromatinModelEnsemble(
            name=models[0].name,
            format=models[0].format,
            first_bin=models[0].first_bin,
            last_bin=models[0].last_bin,
            resolution=models[0].resolution,
            count=len(models),
            metadata_stack=metadata_stack,
            coordinates_stack=coordinates_stack
        )

    @staticmethod
    def from_head(head: ChromatinModelEnsembleHead, coordinates_stack: NDArray[Shape["*, *, 3"], np.float32]) -> 'ChromatinModelEnsemble':
        return ChromatinModelEnsemble(
            name=head.name,
            format=head.format,
            count=head.count,
            first_bin=head.first_bin,
            last_bin=head.last_bin,
            resolution=head.resolution,
            metadata_stack=head.metadata_stack,
            coordinates_stack=coordinates_stack
        )
