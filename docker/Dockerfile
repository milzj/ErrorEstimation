FROM ubuntu:22.04

ENV DEB_PYTHON_INSTALL_LAYOUT=deb_system
ARG DEBIAN_FRONTEND=noninteractive

# python 3.8.11
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y build-essential && \
    apt install -y python3.8 && \
    apt install -y git wget && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

ARG REPO_BRANCH="main"

WORKDIR /repo

# git clone and conda create
RUN git clone --branch ${REPO_BRANCH} --single-branch https://github.com/milzj/ErrorEstimation
RUN cd ErrorEstimation && conda env create -f environment.yml
RUN cd && conda init bash && exec bash





