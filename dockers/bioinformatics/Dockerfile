FROM debian:buster-slim

ENV HAIL_VERSION 0.2.73.devc6f6f09cec08
ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt-get update && apt-get install -y git wget bash bzip2 zip tmux && \
    wget -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 && chmod +x ./jq && cp jq /usr/bin && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c cpg -c bioconda -c conda-forge \
        hail=$HAIL_VERSION \
        cpg-utils \
        google-cloud-sdk \
        dragmap \
        bwa \
        bazam \
        biobambam \
        samtools \
        picard \
        gatk4 \
        bcftools && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c vladsaveliev -c conda-forge tsvfmt && \
    rm -r /root/micromamba/pkgs && \
    wget https://github.com/brentp/somalier/releases/download/v0.2.15/somalier && \
    mv ./somalier /root/micromamba/bin/somalier && \
    chmod +x /root/micromamba/bin/somalier
