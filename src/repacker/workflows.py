from datetime import timedelta

from temporalio import workflow
from temporalio.common import RetryPolicy

from repacker.models import Repack3dgnomeModelEnsembleWorkflowInput, Repack3dgnomeModelEnsembleActivityInput

with workflow.unsafe.imports_passed_through():
    from .activities import repack_3dgnome_model_ensemble_from_bucket


@workflow.defn(name="repack-3dgnome-model-ensemble")
class Repack3dgnomeModelEnsembleWorkflow:

    @workflow.run
    async def run(self, input: Repack3dgnomeModelEnsembleWorkflowInput) -> None:
        await workflow.execute_local_activity(
            activity=repack_3dgnome_model_ensemble_from_bucket,
            arg=Repack3dgnomeModelEnsembleActivityInput(
                source_data_path=input.source_data_path,
                target_data_path=input.target_data_path,
                source_model_name=input.source_model_name,
            ),
            schedule_to_close_timeout=timedelta(minutes=5),
            retry_policy=RetryPolicy(
                maximum_attempts=3,
                initial_interval=timedelta(seconds=1),
                maximum_interval=timedelta(seconds=10),
                non_retryable_error_types=['ValueError'],
            )
        )
