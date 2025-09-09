ARG UBUNTU_VER=22.04

#######################################################################
#     Basic image
#######################################################################
FROM ubuntu:${UBUNTU_VER}

LABEL maintainers="Michele Berselli (berselli.michele@gmail.com)"

## Use latest conda, build for x86_64 and Python 3.9 (default)
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64
ARG PYTHON_VERSION=3.9.17
ARG JDK_VERSION=17.0.3

ENV DEBIAN_FRONTEND noninteractive

#######################################################################
#     General updates & installing necessary Linux components
#######################################################################
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    g++ \
    git \
    less \
    libncurses-dev \
    make \
    time \
    unzip \
    vim \
    wget \
    curl \
    zlib1g-dev \
    liblz4-tool \
    libbz2-dev \
    liblzma-dev \
    locales \
    parallel \
    autoconf \
    autogen \
    pkg-config \
    tabix \
    libcurl4-openssl-dev \
    libssl-dev  

#######################################################################
#     Conda
#######################################################################
## Install Conda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh && \
    bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda3 -b
ENV PATH=/miniconda3/bin:$PATH
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda update -y conda \
    && rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#######################################################################
#     Installing granite and its dependencies
#######################################################################
## Working directory
WORKDIR /usr/local/bin

## Install Python and OpenJDK
RUN conda install -c bioconda -y \
    python=${PYTHON_VERSION} \
    openjdk=${JDK_VERSION}

## Install granite
RUN pip install poetry
RUN poetry config virtualenvs.create false

# Clone and install Granite with poetry
RUN git clone https://github.com/dbmi-bgm/granite.git && \
    cd granite && \
    poetry install --no-root && \
    poetry build && \
    pip install dist/*.whl && \
    cd ..

## Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure && \
    make && \
    cd .. && \
    ln -s samtools-1.17 samtools

## Add mpileup_parser parallel script
COPY mpileupCounts-parallel.sh .
RUN chmod +x mpileupCounts-parallel.sh

#######################################################################
#     Setting up environment
#######################################################################
## Update PATH
ENV PATH=/usr/local/bin/samtools/:$PATH
ENV PATH=/usr/local/bin/granite/:$PATH

## Supporting UTF-8
RUN locale-gen "en_US.UTF-8" && update-locale LC_ALL="en_US.UTF-8"
ENV LC_ALL=en_US.UTF-8

CMD ["bash"]
