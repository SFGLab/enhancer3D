from datetime import timedelta

from temporalio import workflow

from servant.models import Query, Response
from utils.workflow_utils import get_default_retry_policy

with workflow.unsafe.imports_passed_through():
    from .activities import execute_query


@workflow.defn(name="execute-query")
class ExecuteQueryWorkflow:

    @workflow.run
    async def run(self, query: Query) -> Response:
        return await workflow.execute_activity(
            activity=execute_query,
            arg=query,
            start_to_close_timeout=timedelta(hours=6),
            retry_policy=get_default_retry_policy()
        )
