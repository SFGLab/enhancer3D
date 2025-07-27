from typing import Optional, Union

from pydantic import BaseModel, Field

from servant.models import Response, CellLineWithLinksInput, CellLinkCrossComparisonInput


class GetQueryStatusResponse(BaseModel):
    request_id: str
    status: str
    result: Optional[Response] = None


class RequestQuery(BaseModel):
    input: Union[CellLineWithLinksInput, CellLinkCrossComparisonInput] = Field(discriminator="type")


class RequestQueryResponse(BaseModel):
    request_id: str
