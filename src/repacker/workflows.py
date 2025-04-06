import asyncio
from datetime import timedelta

from temporalio import workflow

from repacker.models import Repack3dgnomeModelEnsembleWorkflowInput, Repack3dgnomeModelEnsembleActivityInput, Repack3dgnomeManyModelEnsemblesWorkflowInput, ListAllModelsInBucketActivityInput
from utils.workflow_utils import get_default_retry_policy

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
                source_model_ensemble_name=input.source_model_ensemble_name,
                source_model_name=input.source_model_name,
            ),
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )


@workflow.defn(name="repack-many-3dgnome-model-ensembles")
class Repack3dgnomeManyModelEnsemblesWorkflow:

    @workflow.run
    async def run(self, input: Repack3dgnomeManyModelEnsemblesWorkflowInput) -> None:
        list_all_models_in_bucket_response = await workflow.execute_local_activity(
            activity=list_all_models_in_bucket,
            arg=ListAllModelsInBucketActivityInput(base_paths=input.source_paths),
            start_to_close_timeout=timedelta(hours=2),
            retry_policy=get_default_retry_policy()
        )

        repacking_activities = [
            workflow.execute_local_activity(
                activity=repack_3dgnome_model_ensemble_from_bucket,
                arg=Repack3dgnomeModelEnsembleActivityInput(
                    source_data_path=model_path,
                    target_data_path=input.target_data_path,
                    source_model_ensemble_name="_".join(filter(None, model_path.split("/"))),
                    source_model_name=None,
                ),
                start_to_close_timeout=timedelta(minutes=30),
                retry_policy=get_default_retry_policy()
            )
            for model_path in list_all_models_in_bucket_response.model_paths
        ]

        await asyncio.gather(*repacking_activities)
