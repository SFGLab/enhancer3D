import asyncio
import logging

import dotenv

from calculator.activities import find_potential_pairs_of_enhancers_promoters_for_project, calculate_distances_for_enhancer_promoters_chunk, upsert_project_configuration, persist_distances_for_enhancer_promoters_chunk
from calculator.workflows import CalculateDistancesForProjectWorkflow
from utils.temporal_utils import get_temporal_worker

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


async def main() -> None:
    temporal_worker = await get_temporal_worker(
        task_queue="calculator-task-queue",
        activities=[
            upsert_project_configuration,
            find_potential_pairs_of_enhancers_promoters_for_project,
            calculate_distances_for_enhancer_promoters_chunk,
            persist_distances_for_enhancer_promoters_chunk
        ],
        workflows=[
            CalculateDistancesForProjectWorkflow
        ],
    )

    await temporal_worker.run()


if __name__ == '__main__':
    dotenv.load_dotenv()
    logging.basicConfig()
    asyncio.run(main())
