FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

COPY requirements.txt requirements-dev.txt ./
RUN pip install -r requirements.txt

COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
COPY cpg_workflows README.md setup.py ./
RUN pip install .
