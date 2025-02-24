FROM python:3.10-slim
# xiang2019/spliceq:v1.0.0

# Install Python packages
COPY . /app
WORKDIR /app
RUN pip install -e .

RUN SPLICE-q.py --help