FROM condaforge/mambaforge:4.11.0-0

COPY ["environment.yml", "./"]

RUN apt update && \
    DEBIAN_FRONTEND=noninteractive apt install -y --no-install-recommends tzdata procps nano && \
    apt-get clean

RUN mamba env update --name base --file environment.yml