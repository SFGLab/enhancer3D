from typing import Optional

from pydantic import BaseModel, Field

from servant.models import Response, QueryInput


class GetQueryStatusResponse(BaseModel):
    request_id: str
    status: str
    result: Optional[Response] = None


class RequestQuery(BaseModel):
    input: QueryInput = Field(discriminator="type")


class RequestQueryResponse(BaseModel):
    request_id: str
