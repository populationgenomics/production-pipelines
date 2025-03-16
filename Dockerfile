FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

# installs zstd, used in massively speeding up directory compression
RUN apt update && \
    apt install -y --no-install-recommends \
    zstd && \
    rm -r /var/lib/apt/lists/* && \
    rm -r /var/cache/apt/*

RUN pip install metamist
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
