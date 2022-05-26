FROM debian:buster-slim

ENV VERSION 1.3.0

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt-get update && apt-get install -y \
    git wget bash bzip2 zip \
    build-essential libboost-all-dev \
    && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge \
    bazam biobambam samtools google-cloud-sdk \
    && \
    rm -r /root/micromamba/pkgs

RUN apt update && apt install -y build-essential curl libboost-all-dev && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    git clone --depth 1 --branch $VERSION https://github.com/Illumina/DRAGMAP && \
    cd DRAGMAP && \
    HAS_GTEST=0 make -j 8 && \
    mkdir ${MAMBA_ROOT_PREFIX}/bin && \
    mv build/release/dragen-os ${MAMBA_ROOT_PREFIX}/bin/ && \
    cd .. && \
    rm -r /DRAGMAP
