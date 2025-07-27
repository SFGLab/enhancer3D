from typing import Literal, Union

from pydantic import BaseModel, Field


class CellLineWithLinksInput(BaseModel):
    type: Literal["cell_line_with_links"] = "cell_line_with_links"
    cell_line: str


class CellLinkCrossComparisonInput(BaseModel):
    type: Literal["cell_link_cross_comparison"] = "cell_link_cross_comparison"
    cell_line_base: str
    cell_line_compare: str


class Query(BaseModel):
    request_id: str
    input: Union[CellLineWithLinksInput, CellLinkCrossComparisonInput] = Field(discriminator="type")


class Response(BaseModel):
    request_id: str
    path: str
