FROM python:3.12-slim

WORKDIR /app

# Install requirements
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the source code
COPY src/ src/

# Run the application
ENV PYTHONUNBUFFERED=1
CMD ["python", "src/app_repacker.py"]
