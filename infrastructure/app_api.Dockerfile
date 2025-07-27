FROM python:3.12-slim

EXPOSE 80
WORKDIR /app

# Install build environment
RUN apt-get update -qq && apt-get install -y build-essential gcc zlib1g-dev bedtools

# Install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir fastapi[standard]

# Copy the source code
COPY src/ src/

# Run the application
ENV PYTHONUNBUFFERED=1
CMD ["fastapi", "run", "src/app_api.py", "--port", "80"]
