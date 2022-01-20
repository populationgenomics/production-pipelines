FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:8299b9cf275f677694bb934872ba26f871110d61-hail-0.2.73.devc6f6f09cec08

RUN apt-get update && apt-get install -y tmux

RUN micromamba install -y --prefix $MAMBA_ROOT_PREFIX -c vladsaveliev -c bioconda -c conda-forge \
    google-cloud-sdk \
    dragmap \
    bwa \
    bazam \
    biobambam \
    samtools \
    picard \
    gatk4 \
    bcftools \
    tsvfmt \
    && \
    rm -r /root/micromamba/pkgs
