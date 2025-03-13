import os
from functools import cache
from typing import Optional

from pymongo import MongoClient


@cache
def get_mongo_client(
    connection_string: Optional[str] = None,
) -> MongoClient:
    if not connection_string:
        connection_string = os.getenv("MONGO_CONNECTION_STRING", "mongodb://mongo:27017/")

    return MongoClient(connection_string)
