import asyncio
import logging
from concurrent.futures import ThreadPoolExecutor

import dotenv

from calculator.activities import find_potential_pairs_of_enhancers_promoters_for_project, calculate_distances_for_enhancer_promoters_chunk, upsert_project_configuration
from calculator.workflows import CalculateDistancesForProjectWorkflow
from utils.temporal_utils import get_temporal_worker

logger = logging.getLogger(__name__)


async def main() -> None:
    with ThreadPoolExecutor() as executor:
        temporal_worker = await get_temporal_worker(
            executor=executor,
            task_queue="calculator-task-queue",
            activities=[
                upsert_project_configuration,
                find_potential_pairs_of_enhancers_promoters_for_project,
                calculate_distances_for_enhancer_promoters_chunk,
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
