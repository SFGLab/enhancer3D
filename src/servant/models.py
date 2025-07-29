from typing import Literal, Union, Optional, List

from pydantic import BaseModel, Field


class PartialChromatinRegion(BaseModel):
    chromosome: str
    start: Optional[int] = None
    end: Optional[int] = None


class CellLineWithLinksInput(BaseModel):
    type: Literal["cell_line_with_links"] = "cell_line_with_links"
    cell_line: str
    regions: List[PartialChromatinRegion] = []
    gene_ids: List[str] = []


class CellLinkCrossComparisonInput(BaseModel):
    type: Literal["cell_link_cross_comparison"] = "cell_link_cross_comparison"
    cell_line_base: str
    cell_line_compare: str
    regions: List[PartialChromatinRegion] = []
    gene_ids: List[str] = []


QueryInput = Union[CellLineWithLinksInput, CellLinkCrossComparisonInput]


class Query(BaseModel):
    request_id: str
    input: QueryInput = Field(discriminator="type")


class Response(BaseModel):
    request_id: str
    path: str
