FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:a15ed15c0ca0b7304e0ceca1b873d734720c5d7a-hail-abaf45c99b32250dfc55a1b26a6a0635760e0a99
# this pinned image is a temporary workaround for production-pipelines/issues/477

RUN pip install metamist 
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
