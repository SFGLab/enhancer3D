from datetime import timedelta

from temporalio import workflow

from utils.workflow_utils import get_default_retry_policy, execute_activities_in_batches

with workflow.unsafe.imports_passed_through():
    from .models import CalculateDistancesForProjectWorkflowInput, CalculateDistancesForProjectWorkflowOutput, CalculateDistancesForEnhancerPromotersChunkActivityInput, FindPotentialPairsOfEnhancersPromotersForProjectActivityInput, UpsertProjectConfigurationActivityInput, PersistDistancesForEnhancerPromotersChunkActivityInput
    from .activities import find_potential_pairs_of_enhancers_promoters_for_project, calculate_distances_for_enhancer_promoters_chunk, upsert_project_configuration, persist_distances_for_enhancer_promoters_chunk


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

        find_activity_outputs = await execute_activities_in_batches(
            activity=find_potential_pairs_of_enhancers_promoters_for_project,
            activity_args_list=[
                FindPotentialPairsOfEnhancersPromotersForProjectActivityInput(
                    project=input.project,
                    dataset=dataset,
                    configuration=input.configuration
                )
                for dataset in input.datasets
            ],
            batch_size=100,
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        calculation_outputs = await execute_activities_in_batches(
            activity=calculate_distances_for_enhancer_promoters_chunk,
            activity_args_list=[
                CalculateDistancesForEnhancerPromotersChunkActivityInput(
                    project=input.project,
                    dataset=output.dataset,
                    enhancers_promoters_chunk_path=chunk_path
                )
                for output in find_activity_outputs
                for chunk_path in output.enhancers_promoters_chunk_paths
            ],
            batch_size=100,
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        await execute_activities_in_batches(
            activity=persist_distances_for_enhancer_promoters_chunk,
            activity_args_list=[
                PersistDistancesForEnhancerPromotersChunkActivityInput(
                    project=input.project,
                    dataset=output.dataset,
                    distances_chunk_path=output.distances_chunk_path
                )
                for output in calculation_outputs
            ],
            batch_size=100,
            start_to_close_timeout=timedelta(minutes=30),
            retry_policy=get_default_retry_policy()
        )

        distances_chunk_paths = [output.distances_chunk_path for output in calculation_outputs]
        return CalculateDistancesForProjectWorkflowOutput(distances_chunk_paths=distances_chunk_paths)
