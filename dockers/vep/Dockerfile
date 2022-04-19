FROM debian:stretch-slim

ENV MAMBA_ROOT_PREFIX /root/micromamba
ENV PATH $MAMBA_ROOT_PREFIX/bin:$PATH

RUN apt-get update && \
    apt-get install -y git wget bash bzip2 zip lsb-release curl gnupg2 && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/* && \
    wget -qO- https://api.anaconda.org/download/conda-forge/micromamba/0.8.2/linux-64/micromamba-0.8.2-he9b6cbd_0.tar.bz2 | tar -xvj -C /usr/local bin/micromamba && \
    mkdir $MAMBA_ROOT_PREFIX && \
    micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c bioconda -c conda-forge \
    ensembl-vep=105.0 \
    # Loftee deps:
    samtools perl-list-moreutils perl-dbd-sqlite \
    google-cloud-sdk && \
    rm -r /root/micromamba/pkgs

# Install Bio::Perl for Loftee (conda's perl-bioperl doesn't work)
RUN vep_install --AUTO a --NO_HTSLIB

# Install gcsfuse to mount VEP cache fiels
RUN export GCSFUSE_REPO=gcsfuse-$(lsb_release -c -s) && \
    echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y gcsfuse
