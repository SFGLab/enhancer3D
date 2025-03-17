FROM bitnami/spark:3.5.3 AS base

# Install Livy
ARG LIVY_VERSION=0.8.0-incubating
ARG SCALA_VERSION=2.12
ARG LIVY_HOME=/opt/livy

USER root
RUN install_packages unzip curl
RUN curl "https://downloads.apache.org/incubator/livy/${LIVY_VERSION}/apache-livy-${LIVY_VERSION}_${SCALA_VERSION}-bin.zip" -o "apache-livy-${LIVY_VERSION}-bin.zip"
RUN unzip "apache-livy-${LIVY_VERSION}-bin.zip" \
    && rm -rf "apache-livy-${LIVY_VERSION}-bin.zip" \
    && mv "apache-livy-${LIVY_VERSION}_${SCALA_VERSION}-bin" $LIVY_HOME \
    && mkdir $LIVY_HOME/logs \
    && chown -R 1001:1001 $LIVY_HOME

USER 1001
WORKDIR $LIVY_HOME

COPY --chown=1001:1001 infrastructure/config/livy/livy.conf conf/livy.conf
COPY --chown=1001:1001 infrastructure/config/livy/log4j.properties conf/log4j.properties

# Expose Livy port
EXPOSE 8998

# Force java to open all modules
ENV _JAVA_OPTIONS --add-opens=java.base/java.util=ALL-UNNAMED\
 --add-opens=java.base/java.util.concurrent=ALL-UNNAMED\
 --add-opens=java.base/java.util.concurrent.atomic=ALL-UNNAMED\
 --add-opens=java.base/java.lang=ALL-UNNAMED\
 --add-opens=java.base/java.lang.invoke=ALL-UNNAMED\
 --add-opens=java.base/java.net=ALL-UNNAMED\
 --add-opens=java.base/sun.nio.ch=ALL-UNNAMED\
 --add-opens=java.base/jdk.internal.loader=ALL-UNNAMED\
 --add-opens=java.base/java.security=ALL-UNNAMED

# Start Livy
CMD ["sh", "-c", "bin/livy-server", "start"]
