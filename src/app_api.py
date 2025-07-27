import uuid
from typing import Optional

from fastapi import FastAPI

from api.models import GetQueryStatusResponse, RequestQuery, RequestQueryResponse
from utils.temporal_utils import get_temporal_client

app = FastAPI()

@app.get("/query/{id}")
async def get_query_status(id: str) -> Optional[GetQueryStatusResponse]:
    temporal_client = await get_temporal_client()
    workflow_run_handle = temporal_client.get_workflow_handle(
        workflow_id="execute-query",
        run_id=id
    )

    workflow_run_description = await workflow_run_handle.describe()
    workflow_status = workflow_run_description.status.name if workflow_run_description.status else "UNKNOWN"

    if workflow_status != "COMPLETED":
        return GetQueryStatusResponse(
            request_id=workflow_run_description.run_id,
            status=workflow_status
        )

    workflow_result = await workflow_run_handle.result()
    return GetQueryStatusResponse(
        request_id=workflow_run_description.run_id,
        status=workflow_status,
        result=workflow_result
    )


@app.post("/query")
async def request_query(request: RequestQuery) -> RequestQueryResponse:
    temporal_client = await get_temporal_client()

    workflow_run_id = str(uuid.uuid4())
    await temporal_client.start_workflow(
        workflow="execute-query",
        arg=request.input,
        id=workflow_run_id,
        task_queue="servant-task-queue",
    )

    return RequestQueryResponse(request_id=workflow_run_id)
