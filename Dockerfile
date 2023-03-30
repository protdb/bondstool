FROM ubuntu:22.04
USER root
WORKDIR /app
RUN apt-get update
RUN apt-get install -y python3-pip python-is-python3
RUN apt-get install -y pymol
COPY ./requirements.txt /app/requirements.txt
RUN pip install -r /app/requirements.txt
ADD . /app
CMD ["python", "-u", "run_structures_set.py"]