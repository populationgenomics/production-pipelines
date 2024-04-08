FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
COPY README.md setup.py cpg_workflows .
RUN pip install .
