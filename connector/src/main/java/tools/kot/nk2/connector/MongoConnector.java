package tools.kot.nk2.connector;

import com.fasterxml.jackson.databind.ObjectMapper;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoClients;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.api.java.UDF2;
import org.apache.spark.sql.catalyst.expressions.GenericRow;
import org.apache.spark.sql.catalyst.expressions.GenericRowWithSchema;
import org.apache.spark.sql.types.*;
import org.bson.Document;
import scala.collection.Seq;
import scala.jdk.CollectionConverters;
import java.io.Serial;
import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.StreamSupport;

public class MongoConnector implements Serializable, UDF2<String, Seq<GenericRowWithSchema>, Seq<GenericRowWithSchema>> {

    @Serial
    private static final long serialVersionUID = 1L;

    private static final ObjectMapper OBJECT_MAPPER = new ObjectMapper();

    @Override
    public Seq<GenericRowWithSchema> call(String configurationJson, Seq<GenericRowWithSchema> conditions) throws Exception {
        final var config = OBJECT_MAPPER.readValue(configurationJson, MongoConnectorConfiguration.class);

        if(!(config.structure() instanceof StructType structType)) {
            throw new IllegalArgumentException("Configuration structure must be a StructType, instead got: " + config.structure());
        }

        // Connect to MongoDB
        try (MongoClient client = MongoClients.create(config.connectionUri())) {
            final var database = client.getDatabase(config.databaseName());
            final var collection = database.getCollection(config.collectionName());

            final var conditionsIterable = CollectionConverters.asJavaIterable(conditions);

            // Convert conditions to MongoDB query format
            final var criteriaDocuments = new ArrayList<Document>();
            for (final var row : conditionsIterable) {
                final var criteriaDoc = new Document();
                for (final var field : row.schema().fields()) {
                    criteriaDoc.append(field.name(), row.getAs(field.name()));
                }

                criteriaDocuments.add(criteriaDoc);
            }

            // Create the MongoDB query with the specified root condition
            final var resolvedRootCondition = config.rootCondition() == null || config.rootCondition().isEmpty()
                ? "or" // Default to "or" if no root condition is specified
                : config.rootCondition();

            final var query = new Document("$" + resolvedRootCondition, criteriaDocuments);

            // Execute the query with or without projection
            Document projectionDocument = null;
            if (config.projection() != null && !config.projection().isEmpty()) {
                projectionDocument = new Document();
                for (final var field : config.projection()) {
                    projectionDocument.append(field, 1);
                }
            }

            final var results = StreamSupport.stream(collection.find(query).projection(projectionDocument).spliterator(), true)
                .map(document -> coerceDocumentToSchema(document, structType))
                .toList();

            return CollectionConverters.asScalaBuffer(results);
        }
    }

    private static GenericRowWithSchema coerceDocumentToSchema(Document document, StructType schema) {
        final var values = new Object[schema.fields().length];
        for (int i = 0; i < schema.fields().length; i++) {
            final var field = schema.fields()[i];
            final var value = document.get(field.name());

            final var dataType = field.dataType();
            if (dataType instanceof StructType structDataType) {
                final var documentValue = (Document) value;
                values[i] = coerceDocumentToSchema(documentValue, structDataType);
            } else if (dataType instanceof ArrayType arrayDataType) {
                final var listValue = (List<?>) value;
                values[i] = coerceArrayToSchema(listValue, arrayDataType);
            } else {
                values[i] = coercePrimitiveToSchema(value, dataType);
            }
        }

        return new GenericRowWithSchema(values, schema);
    }

    private static Object[] coerceArrayToSchema(List<?> list, ArrayType arrayType) {
        final var elementType = arrayType.elementType();
        final var values = new Object[list.size()];

        for (int i = 0; i < list.size(); i++) {
            final var item = list.get(i);
            if (item instanceof Document document) {
                values[i] = coerceDocumentToSchema(document, (StructType) elementType);
            } else if (item instanceof List<?> itemList) {
                values[i] = coerceArrayToSchema(itemList, (ArrayType) elementType);
            } else {
                values[i] = coercePrimitiveToSchema(item, elementType);
            }
        }

        return values;
    }

    private static Object coercePrimitiveToSchema(Object primitive, DataType dataType) {
        if (dataType instanceof StringType) {
            return primitive.toString();
        } else if (dataType instanceof IntegerType) {
            return Integer.parseInt(primitive.toString());
        } else if (dataType instanceof LongType) {
            return Long.parseLong(primitive.toString());
        } else if (dataType instanceof DoubleType) {
            return Double.parseDouble(primitive.toString());
        } else if (dataType instanceof FloatType) {
            return Float.parseFloat(primitive.toString());
        } else if (dataType instanceof ShortType) {
            return Short.parseShort(primitive.toString());
        } else if (dataType instanceof ByteType) {
            return Byte.parseByte(primitive.toString());
        } else if (dataType instanceof DecimalType decimalType) {
            return new BigDecimal(primitive.toString()).setScale(decimalType.scale(), java.math.RoundingMode.HALF_UP);
        } else if (dataType instanceof BooleanType) {
            return Boolean.parseBoolean(primitive.toString());
        } else {
            throw new IllegalArgumentException("Unsupported primitive type: " + dataType);
        }
    }
}
