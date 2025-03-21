FROM bitnami/spark:3.5.3 AS base

USER root

# Install build environment
RUN apt-get update -qq && apt-get install -y build-essential gcc zlib1g-dev bedtools

# Install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
