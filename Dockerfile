FROM python:3.11-slim

WORKDIR /project
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .

CMD ["python", "experiments/run_smoke_experiment.py"]
