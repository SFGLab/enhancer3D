from typing import Sequence

from pymongo import UpdateOne

from utils.pydantic_utils import BaseDatabaseModel
from utils.mongo_utils import get_mongo_client


def upsert_many_database_models(
    database_name: str,
    collection_name: str,
    data: Sequence[BaseDatabaseModel],
    batch_size: int = 500
) -> None:
    mongo_client = get_mongo_client()
    database = mongo_client.get_database(database_name)
    collection = database.get_collection(collection_name)

    for batch_position in range(0, len(data), batch_size):
        batch_id = batch_position // batch_size
        data_batch = data[batch_position:batch_position + batch_size]

        collection.bulk_write([
            UpdateOne(
                filter={
                    "_id": entry.model_dump(include=set(entry.index))
                },
                update={
                    "$set": {
                        "_id": entry.model_dump(include=set(entry.index)),
                        **entry.model_dump(include=set(entry.fields))
                    }
                },
                upsert=True,
            )
            for entry in data_batch
        ])
