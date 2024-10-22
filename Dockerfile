
FROM python:3.10-slim-buster

WORKDIR /python-docker
RUN apt-get update \
    && apt-get -y install libpq-dev gcc \
    && pip install psycopg2\
    && pip install minio
    

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .
EXPOSE 3000
CMD [ "python", "app.py"]