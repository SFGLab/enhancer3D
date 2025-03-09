import asyncio
import logging
from concurrent.futures import ThreadPoolExecutor

import dotenv

from repacker.activities import repack_3dgnome_model_ensemble_from_bucket
from repacker.workflows import Repack3dgnomeModelEnsembleWorkflow
from utils.temporal_utils import get_temporal_worker

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


async def main() -> None:
    with ThreadPoolExecutor() as executor:
        temporal_worker = await get_temporal_worker(
            executor=executor,
            task_queue="repacker-task-queue",
            activities=[repack_3dgnome_model_ensemble_from_bucket],
            workflows=[Repack3dgnomeModelEnsembleWorkflow],
        )

        await temporal_worker.run()


if __name__ == '__main__':
    dotenv.load_dotenv()
    asyncio.run(main())
