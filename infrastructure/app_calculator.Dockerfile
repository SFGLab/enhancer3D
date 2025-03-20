FROM python:3.12-slim

WORKDIR /app

# Install build environment
RUN apt-get update -qq && apt-get install -y build-essential gcc zlib1g-dev

# Install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the source code
COPY src/ src/

# Run the application
ENV PYTHONUNBUFFERED=1
CMD ["python", "src/app_calculator.py"]
