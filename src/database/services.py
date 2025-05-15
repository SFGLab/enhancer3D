from typing import Sequence

from pydantic.fields import FieldInfo
from pymongo import UpdateOne

from utils.mongo_utils import get_mongo_client
from utils.pydantic_utils import BaseDatabaseModel
from utils.scylla_utils import get_scylla_session


def upsert_many_database_models_mongo(
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


def upset_many_database_models_cassandra(
    keyspace: str,
    table: str,
    data: Sequence[BaseDatabaseModel],
    batch_size: int = 500
) -> None:
    session = get_scylla_session(keyspace=keyspace)

    # Create the table if it doesn't exist with the correct schema
    def type_by_field_info(model_field_info: FieldInfo):
        if model_field_info.annotation == str:
            return 'text'
        elif model_field_info.annotation == int:
            return 'int'
        elif model_field_info.annotation == float:
            return 'float'
        elif model_field_info.annotation == bool:
            return 'boolean'
        elif model_field_info.annotation == list:
            return f'list<{type(model_field_info.annotation.__type_params__[0]).__name__}>'
        else:
            raise ValueError(f"Unsupported type: {model_field_info.annotation}")

    columns = ', '.join(
        f"{field.name} {type_by_field_info(field)}"
        for field in data[0].model_fields.values()
    )
    index_columns = ', '.join(data[0].index)
    query = f"""
        CREATE TABLE IF NOT EXISTS {table} (
            {columns},
            PRIMARY KEY ({index_columns})
        )
    """
    session.execute(query)

    # Insert data in batches
    for batch_position in range(0, len(data), batch_size):
        batch_id = batch_position // batch_size
        data_batch = data[batch_position:batch_position + batch_size]

        for entry in data_batch:
            query = f"""
                INSERT INTO {table} ({', '.join(entry.fields)})
                VALUES ({', '.join(['%s'] * len(entry.fields))})
            """
            values = [entry.model_dump(include=set(entry.fields))]
            session.execute(query, values)
