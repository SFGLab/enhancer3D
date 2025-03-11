import os
from typing import Optional

from pymongo import MongoClient


def get_mongo_client(
    connection_string: Optional[str] = None,
) -> MongoClient:
    if connection_string:
        connection_string = os.getenv("MONGO_CONNECTION_STRING", "mongodb://mongo:27017/")

    return MongoClient(connection_string)
