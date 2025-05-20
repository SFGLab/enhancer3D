from datetime import datetime
from typing import Sequence, get_args, Type

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
                    "_id": entry.model_dump(include=set(entry.index_fields))
                },
                update={
                    "$set": {
                        "_id": entry.model_dump(include=set(entry.index_fields)),
                        **entry.model_dump(include=set(entry.non_index_fields))
                    }
                },
                upsert=True,
            )
            for entry in data_batch
        ])


def upsert_many_database_models_cassandra(
    keyspace: str,
    table: str,
    data: Sequence[BaseDatabaseModel],
    batch_size: int = 500
) -> None:
    session = get_scylla_session(keyspace=keyspace)

    # Create the table if it doesn't exist with the correct schema
    def cassandra_type(python_type: Type):
        if python_type == str:
            return 'text'
        elif python_type == int:
            return 'int'
        elif python_type == float:
            return 'float'
        elif python_type == bool:
            return 'boolean'
        elif python_type.__name__ == 'List':
            type_args = get_args(python_type)
            return f'list<{cassandra_type(type_args[0])}>'
        elif python_type == datetime:
            return 'timestamp'
        else:
            raise ValueError(f"Unsupported type: {python_type}")

    def create_cassandra_table(example: BaseDatabaseModel):
        columns = ', '.join(
            f"{name} {cassandra_type(field.annotation)}"
            for name, field in example.model_fields.items()
        )

        index_columns = ', '.join(example.index_fields)
        query = f"""
            CREATE TABLE IF NOT EXISTS {table} (
                {columns},
                PRIMARY KEY ({index_columns})
            )
        """

        session.execute(query)

    def insert_cassandra_data(entry: BaseDatabaseModel):
        query = f"""
            INSERT INTO {table} ({', '.join(entry.all_fields)})
            VALUES ({', '.join(['%s'] * len(entry.all_fields))})
        """

        model_dict = entry.model_dump(include=set(entry.all_fields))
        values = [model_dict.get(field) for field in entry.all_fields]
        session.execute(query, values)

    # Create the table with the first entry as an example
    create_cassandra_table(data[0])

    # Insert data in batches
    for batch_position in range(0, len(data), batch_size):
        batch_id = batch_position // batch_size
        data_batch = data[batch_position:batch_position + batch_size]

        for entry in data_batch:
            insert_cassandra_data(entry)
