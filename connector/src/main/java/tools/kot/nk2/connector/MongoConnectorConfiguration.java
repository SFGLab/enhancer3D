package tools.kot.nk2.connector;

import org.apache.spark.sql.types.DataType;

import java.util.List;

public record MongoConnectorConfiguration(
    String connectionUri,
    String databaseName,
    String collectionName,
    String rootCondition,
    List<String> projection,
    DataType structure
) {}
