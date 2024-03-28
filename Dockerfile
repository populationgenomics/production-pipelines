FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

RUN pip install metamist
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
