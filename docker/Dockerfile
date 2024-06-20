FROM ubuntu:22.04

ENV DEB_PYTHON_INSTALL_LAYOUT=deb_system
ARG DEBIAN_FRONTEND=noninteractive

# python 3.8.11
RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt-get update && \
    apt install -y python3.8 && \
    apt install -y wget git curl nano && \
    rm -rf /var/lib/apt/lists/*

# Miniconda

ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
ARG PY_VER=3.8.11
ARG TF_VER=2.5.0
ARG CONDA_VER
ARG OS_TYPE

# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

ARG REPO_BRANCH="main"

WORKDIR /repo

# Copy pyproject.toml first so that we done need to reinstall in case another file
# is changing after rebuiding docker image
RUN git clone --branch ${REPO_BRANCH} --single-branch https://github.com/milzj/ErrorEstimation
RUN cd ErrorEstimation && conda env create -f environment.yml
RUN cd && conda init bash && exec bash
#RUN conda activate ErrorEstimation
#RUN cd ErrorEstimation && pip install pytest && pytest




