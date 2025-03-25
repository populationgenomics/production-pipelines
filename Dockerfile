FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

ADD requirements.txt .
RUN pip install --no-cache-dir -r  requirements.txt

COPY README.md setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install --no-cache-dir .
