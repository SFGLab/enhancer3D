FROM bitnami/spark:3.5.3 AS base

# Install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
