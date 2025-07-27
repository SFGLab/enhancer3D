from typing import Literal

from pydantic import BaseModel, Field


class Query(BaseModel):
    request_id: str
    input: "CellLineWithLinksInput" = Field(discriminator="type")


class CellLineWithLinksInput(BaseModel):
    type: Literal["cell_line_with_links"] = "cell_line_with_links"
    cell_line: str


class Response(BaseModel):
    request_id: str
    path: str
