from typing import NamedTuple, Dict, Any, List

import numpy as np


ChromatinModelMetadata = Dict[str, Any]


class ChromatinModel(NamedTuple):
    id: int
    name: str
    format: str
    metadata: ChromatinModelMetadata
    coordinates: np.ndarray[np.float64]


class ChromatinModelEnsemble(NamedTuple):
    name: str
    format: str
    count: int

    metadata_stack: List[ChromatinModelMetadata]
    coordinates_stack: np.ndarray[np.float64]

    def get_model(self, model_id: int) -> ChromatinModel:
        if not 0 <= model_id < self.count:
            raise ValueError(f'Model ensemble does not have {model_id} model.')

        coordinates = self.coordinates_stack[model_id]
        metadata = self.metadata_stack[model_id]
        return ChromatinModel(
            id=model_id,
            name=self.name,
            format=self.format,
            metadata=metadata,
            coordinates=coordinates
        )
