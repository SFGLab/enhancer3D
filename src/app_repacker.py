import asyncio
import logging
from concurrent.futures import ThreadPoolExecutor

import dotenv

from repacker.activities import repack_3dgnome_model_ensemble_from_bucket, list_all_models_in_bucket
from repacker.workflows import Repack3dgnomeModelEnsembleWorkflow, Repack3dgnomeManyModelEnsemblesWorkflow
from utils.temporal_utils import get_temporal_worker

logger = logging.getLogger(__name__)


async def main() -> None:
    temporal_worker = await get_temporal_worker(
        task_queue="repacker-task-queue",
        activities=[
            repack_3dgnome_model_ensemble_from_bucket,
            list_all_models_in_bucket,
        ],
        workflows=[
            Repack3dgnomeModelEnsembleWorkflow,
            Repack3dgnomeManyModelEnsemblesWorkflow,
        ],
    )

    await temporal_worker.run()


if __name__ == '__main__':
    dotenv.load_dotenv()
    logging.basicConfig()
    asyncio.run(main())
