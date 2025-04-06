import asyncio
from datetime import timedelta

from temporalio import workflow

from utils.workflow_utils import get_default_retry_policy

with workflow.unsafe.imports_passed_through():
    from .models import CalculateDistancesForProjectWorkflowInput, CalculateDistancesForProjectWorkflowOutput, CalculateDistancesForEnhancerPromotersChunkActivityInput, FindPotentialPairsOfEnhancersPromotersForProjectActivityInput, UpsertProjectConfigurationActivityInput, PersistDistancesForEnhancerPromotersChunkActivityInput, PreloadDatasetsForProjectActivityInput
    from .activities import find_potential_pairs_of_enhancers_promoters_for_project, calculate_distances_for_enhancer_promoters_chunk, upsert_project_configuration, persist_distances_for_enhancer_promoters_chunk, preload_datasets_for_project


@workflow.defn(name="calculate-distances-for-project")
class CalculateDistancesForProjectWorkflow:

    @workflow.run
    async def run(self, input: CalculateDistancesForProjectWorkflowInput) -> CalculateDistancesForProjectWorkflowOutput:
        await workflow.execute_activity(
            activity=upsert_project_configuration,
            arg=UpsertProjectConfigurationActivityInput(
                project=input.project,
                datasets=input.datasets,
                configuration=input.configuration
            ),
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        await workflow.execute_activity(
            activity=preload_datasets_for_project,
            arg=PreloadDatasetsForProjectActivityInput(
                project=input.project,
                datasets=input.datasets,
                configuration=input.configuration
            ),
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        find_activities = [
            workflow.execute_activity(
                activity=find_potential_pairs_of_enhancers_promoters_for_project,
                arg=FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(
                    project=input.project,
                    dataset=dataset,
                    configuration=input.configuration
                ),
                start_to_close_timeout=timedelta(minutes=30),
                retry_policy=get_default_retry_policy()
            )
            for dataset in input.datasets
        ]

        find_activity_outputs = await asyncio.gather(*find_activities)

        calculation_activities = [
            workflow.execute_activity(
                activity=calculate_distances_for_enhancer_promoters_chunk,
                arg=CalculateDistancesForEnhancerPromotersChunkActivityInput(
                    project=input.project,
                    dataset=output.dataset,
                    enhancers_promoters_chunk_path=chunk_path
                ),
                schedule_to_close_timeout=timedelta(minutes=30),
                retry_policy=get_default_retry_policy()
            )
            for output in find_activity_outputs
            for chunk_path in output.enhancers_promoters_chunk_paths
        ]

        calculation_outputs = await asyncio.gather(*calculation_activities)
        persist_activities = [
            workflow.execute_activity(
                activity=persist_distances_for_enhancer_promoters_chunk,
                arg=PersistDistancesForEnhancerPromotersChunkActivityInput(
                    project=input.project,
                    dataset=output.dataset,
                    distances_chunk_path=output.distances_chunk_path
                ),
                start_to_close_timeout=timedelta(minutes=30),
                retry_policy=get_default_retry_policy()
            )
            for output in calculation_outputs
        ]

        distances_chunk_paths = [output.distances_chunk_path for output in calculation_outputs]

        await asyncio.gather(*persist_activities)
        return CalculateDistancesForProjectWorkflowOutput(distances_chunk_paths=distances_chunk_paths)
