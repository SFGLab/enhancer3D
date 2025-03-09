import asyncio
from datetime import timedelta

from temporalio import workflow
from temporalio.common import RetryPolicy

from repacker.models import Repack3dgnomeModelEnsembleWorkflowInput, Repack3dgnomeModelEnsembleActivityInput, Repack3dgnomeManyModelEnsemblesWorkflowInput, ListAllModelsInBucketActivityInput

with workflow.unsafe.imports_passed_through():
    from .activities import repack_3dgnome_model_ensemble_from_bucket, list_all_models_in_bucket


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


@workflow.defn(name="repack-many-3dgnome-model-ensembles")
class Repack3dgnomeManyModelEnsemblesWorkflow:

    @workflow.run
    async def run(self, input: Repack3dgnomeManyModelEnsemblesWorkflowInput) -> None:
        list_all_models_in_bucket_response = await workflow.execute_local_activity(
            activity=list_all_models_in_bucket,
            arg=ListAllModelsInBucketActivityInput(base_paths=input.source_paths),
            schedule_to_close_timeout=timedelta(minutes=5),
            retry_policy=RetryPolicy(
                maximum_attempts=3,
                initial_interval=timedelta(seconds=1),
                maximum_interval=timedelta(seconds=10),
                non_retryable_error_types=['ValueError'],
            )
        )

        repacking_activities = [
            workflow.execute_local_activity(
                activity=repack_3dgnome_model_ensemble_from_bucket,
                arg=Repack3dgnomeModelEnsembleActivityInput(
                    source_data_path=model_path,
                    target_data_path=input.target_data_path,
                    source_model_name=None,
                ),
                schedule_to_close_timeout=timedelta(minutes=5),
                retry_policy=RetryPolicy(
                    maximum_attempts=3,
                    initial_interval=timedelta(seconds=1),
                    maximum_interval=timedelta(seconds=10),
                    non_retryable_error_types=['ValueError'],
                )
            )
            for model_path in list_all_models_in_bucket_response.model_paths
        ]

        await asyncio.gather(*repacking_activities)
