from abc import ABC, abstractmethod
from typing import Sequence, Set

from pydantic import BaseModel


class BaseDatabaseModel(ABC, BaseModel):

    @property
    @abstractmethod
    def index(self) -> Sequence[str]:
        pass

    @property
    def fields(self) -> Set[str]:
        return (
            set(self.model_fields.keys())
            .union(set(self.model_computed_fields.keys()))
            .difference(set(self.index))
        )
