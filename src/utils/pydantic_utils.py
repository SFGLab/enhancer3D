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



class BaseHashableModel(BaseModel):

    def __hash__(self):
        return hash(tuple(self.model_dump().values()))

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return self.model_dump() == other.model_dump()

    def __ne__(self, other):
        return not self.__eq__(other)
