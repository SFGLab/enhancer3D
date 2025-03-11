from datetime import timedelta

from temporalio.common import RetryPolicy


def get_default_retry_policy() -> RetryPolicy:
    return RetryPolicy(
        maximum_attempts=3,
        initial_interval=timedelta(seconds=1),
        maximum_interval=timedelta(seconds=10),
        non_retryable_error_types=['ValueError'],
    )
