from typing import Literal, Union

from pydantic import BaseModel, Field


class CellLineWithLinksInput(BaseModel):
    type: Literal["cell_line_with_links"] = "cell_line_with_links"
    cell_line: str


class Query(BaseModel):
    request_id: str
    input: Union[CellLineWithLinksInput] = Field(discriminator="type")


class Response(BaseModel):
    request_id: str
    path: str
