import os

from pyspark.sql import SparkSession
from temporalio import activity

from servant.models import Query, Response
from servant.queries import query_cell_line_with_links, query_cell_link_cross_comparison

spark = (
    SparkSession.builder
    .appName("servant")
    .master(os.getenv("SPARK_MASTER", "local[*]"))
    .config("spark.cores.max", int(os.getenv("SPARK_CORES_MAX", "4")))
    .config("spark.executor.instances", int(os.getenv("SPARK_EXECUTOR_INSTANCES", "1")))
    .config("spark.executor.memory", os.getenv("SPARK_EXECUTOR_MEMORY", "4g"))
    .config("spark.driver.memory", os.getenv("SPARK_DRIVER_MEMORY", "1g"))
    .config("spark.driver.maxResultSize", os.getenv("SPARK_DRIVER_MAX_RESULT_SIZE", "1g"))
    .config("spark.authenticate", os.getenv("SPARK_RPC_AUTHENTICATION_ENABLED", "false").lower() == "true")
    .config("spark.authenticate.secret", os.getenv("SPARK_RPC_AUTHENTICATION_SECRET", ""))
    .config("spark.network.crypto.enabled", os.getenv("SPARK_RPC_ENCRYPTION", "false").lower() == "true")
    .config("spark.jars.packages", "org.apache.hadoop:hadoop-aws:3.3.2,com.amazonaws:aws-java-sdk-bundle:1.12.262")
    .config("spark.hadoop.fs.s3a.impl", "org.apache.hadoop.fs.s3a.S3AFileSystem")
    .config("spark.hadoop.fs.s3a.aws.credentials.provider", "org.apache.hadoop.fs.s3a.SimpleAWSCredentialsProvider")
    .config("spark.hadoop.fs.s3a.access.key", os.getenv("BUCKET_ACCESS_KEY", ""))
    .config("spark.hadoop.fs.s3a.secret.key", os.getenv("BUCKET_SECRET_KEY", ""))
    .config("spark.hadoop.fs.s3a.endpoint", os.getenv("BUCKET_ENDPOINT", ""))
    .config("spark.hadoop.fs.s3a.path.style.access", os.getenv("SPARK_S3A_PATH_STYLE_ACCESS", "true").lower() == "true")
    .config("spark.hadoop.fs.s3a.fast.upload", os.getenv("SPARK_S3A_FAST_UPLOAD", "true").lower() == "true")
    .getOrCreate()
)


@activity.defn
def execute_query(query: Query) -> Response:
    """
    Execute a query based on the provided input and return the response.
    """
    if query.input.type == "cell_line_with_links":
        return query_cell_line_with_links(spark, query.request_id, query.input)

    if query.input.type == "cell_link_cross_comparison":
        return query_cell_link_cross_comparison(spark, query.request_id, query.input)

    raise ValueError(f"Unsupported query type: {query.input.type}")
