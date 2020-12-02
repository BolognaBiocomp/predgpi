# Base Image
FROM python:3.8-slim-buster

# Metadata
LABEL base.image="python:3.8-slim-buster"
LABEL version="1.0"
LABEL software="PedGPI"
LABEL software.version="202001"
LABEL description="an open source software tool to predict GPI-anchors in proteins"
LABEL website="https://github.com/BolognaBiocomp/predgpi"
LABEL documentation="https://github.com/BolognaBiocomp/predgpi"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

WORKDIR /usr/src/predgpi

COPY requirements.txt .

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    useradd -m pregpi

USER pregpi

COPY . .

WORKDIR /data/

ENV PREDGPI_HOME='/usr/src/predgpi' PATH=/usr/src/predgpi:$PATH

ENTRYPOINT ["/usr/src/predgpi/predgpi.py"]
