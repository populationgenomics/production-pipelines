FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:db0164a7e4eb0ff73d2d078812213962aa7a5da9-hail-1a26f7a4a3d8496042c504153e81b6d7c91ee8af

RUN pip install metamist
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
