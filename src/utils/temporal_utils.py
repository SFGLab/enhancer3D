import logging
import os
from concurrent.futures import Executor
from functools import cache
from typing import Sequence, Callable, Type, Optional

from temporalio.client import Client
from temporalio.contrib.pydantic import pydantic_data_converter
from temporalio.worker import Worker

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


@cache
async def get_temporal_client(
    temporal_endpoint: Optional[str] = None,
    temporal_namespace: Optional[str] = None,
) -> Client:
    if temporal_endpoint is None:
        temporal_endpoint = os.environ.get("TEMPORAL_ENDPOINT", "temporal:7233")

    if temporal_namespace is None:
        temporal_namespace = os.environ.get("TEMPORAL_NAMESPACE", "default")

    logger.info(f"Connecting to Temporal at {temporal_endpoint} in namespace {temporal_namespace}")
    return await Client.connect(
        namespace=temporal_namespace,
        target_host=temporal_endpoint,
        data_converter=pydantic_data_converter
    )


async def get_temporal_worker(
    executor: Executor,
    task_queue: str,
    activities: Sequence[Callable],
    workflows: Sequence[Type],
    temporal_endpoint: Optional[str] = None,
    temporal_namespace: Optional[str] = None,
) -> Worker:
    client = await get_temporal_client(temporal_endpoint, temporal_namespace)
    logger.info(f"Starting worker for task queue {task_queue}")
    return Worker(
        client,
        task_queue=task_queue,
        workflows=workflows,
        activities=activities,
        activity_executor=executor
    )
