import os
from functools import cache
from typing import Optional

from cassandra.cluster import Cluster, Session, ConsistencyLevel


@cache
def get_scylla_session(contact_points: Optional[str] = None, keyspace: Optional[str] = None) -> Session:
    if not contact_points:
        contact_points = os.getenv("SCYLLA_CONTACT_POINTS", "scylla-master:9042")

    cluster = Cluster(contact_points.split(","))

    # If keyspace is provided, and does not exist, create it
    if keyspace:
        session = cluster.connect()
        try:
            session.execute(f"CREATE KEYSPACE IF NOT EXISTS {keyspace} WITH REPLICATION = {{ 'class': 'SimpleStrategy', 'replication_factor': 1 }}")
        except Exception as e:
            print(f"Error creating keyspace: {e}")
        finally:
            session.shutdown()

    session = cluster.connect(keyspace) if keyspace else cluster.connect()
    session.default_consistency_level = ConsistencyLevel.QUORUM

    return session
