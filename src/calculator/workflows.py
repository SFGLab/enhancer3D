import asyncio
from datetime import timedelta

from temporalio import workflow

from calculator.models import CalculateDistancesForProjectWorkflowInput, CalculateDistancesForProjectWorkflowOutput, CalculateDistancesForEnhancerPromotersChunkActivityInput, FindPotentialPairsOfEnhancersPromotersForProjectActivityInput, UpsertProjectConfigurationActivityInput
from utils.workflow_utils import get_default_retry_policy

with workflow.unsafe.imports_passed_through():
    from .activities import find_potential_pairs_of_enhancers_promoters_for_project, calculate_distances_for_enhancer_promoters_chunk, upsert_project_configuration


@workflow.defn(name="calculate-distances-for-project")
class CalculateDistancesForProjectWorkflow:

    @workflow.run
    async def run(self, input: CalculateDistancesForProjectWorkflowInput) -> CalculateDistancesForProjectWorkflowOutput:
        await workflow.execute_local_activity(
            activity=upsert_project_configuration,
            arg=UpsertProjectConfigurationActivityInput(
                project=input.project,
                configuration=input.configuration
            ),
            schedule_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        find_activity_output = await workflow.execute_local_activity(
            activity=find_potential_pairs_of_enhancers_promoters_for_project,
            arg=FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(
                project=input.project,
                configuration=input.configuration
            ),
            schedule_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        calculation_activities = [
            workflow.execute_local_activity(
                activity=calculate_distances_for_enhancer_promoters_chunk,
                arg=CalculateDistancesForEnhancerPromotersChunkActivityInput(
                    project=input.project,
                    enhancers_promoters_chunk_path=chunk_path
                ),
                schedule_to_close_timeout=timedelta(minutes=30),
                retry_policy=get_default_retry_policy()
            )
            for chunk_path in find_activity_output.enhancers_promoters_chunk_paths
        ]

        calculation_outputs = await asyncio.gather(*calculation_activities)
        distances_chunk_paths = [output.distances_chunk_path for output in calculation_outputs]
        return CalculateDistancesForProjectWorkflowOutput(distances_chunk_paths=distances_chunk_paths)
