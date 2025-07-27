import asyncio
import logging

import dotenv

from servant.activities import execute_query
from servant.workflows import ExecuteQueryWorkflow
from utils.temporal_utils import get_temporal_worker

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


async def main() -> None:
    temporal_worker = await get_temporal_worker(
        task_queue="servant-task-queue",
        activities=[
            execute_query,
        ],
        workflows=[
            ExecuteQueryWorkflow,
        ],
    )

    await temporal_worker.run()


if __name__ == '__main__':
    dotenv.load_dotenv()
    logging.basicConfig()
    asyncio.run(main())
