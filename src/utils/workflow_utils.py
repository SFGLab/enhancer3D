import asyncio
from datetime import timedelta
from typing import Any, List, Optional

from temporalio import workflow
from temporalio.common import RetryPolicy


def get_default_retry_policy() -> RetryPolicy:
    return RetryPolicy(
        maximum_attempts=3,
        initial_interval=timedelta(seconds=1),
        maximum_interval=timedelta(seconds=10),
        non_retryable_error_types=['ValueError'],
    )


async def execute_activities_in_batches(
    activity: Any,
    activity_args_list: List[Any],
    batch_size: int = 100,
    start_to_close_timeout: Optional[timedelta] = None,
    retry_policy: Optional[RetryPolicy] = None,
) -> List[Any]:
    results = []

    for i in range(0, len(activity_args_list), batch_size):
        batch_args = activity_args_list[i:i + batch_size]

        # Create activities only for the current batch
        batch_activities = [
            workflow.execute_activity(
                activity=activity,
                arg=arg,
                start_to_close_timeout=start_to_close_timeout,
                retry_policy=retry_policy
            )
            for arg in batch_args
        ]

        # Execute this batch and wait for completion
        batch_results = await asyncio.gather(*batch_activities)
        results.extend(batch_results)

    return results
