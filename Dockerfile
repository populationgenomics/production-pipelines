FROM australia-southeast1-docker.pkg.dev/analysis-runner/images/driver:latest

RUN pip install metamist \
&& pip install git+https://github.com/Illumina/ica-sdk-python.git \
&& pip install typing-extensions --upgrade
COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
COPY seqr-loading-pipelines/hail_scripts hail_scripts
COPY gnomad_methods/gnomad gnomad
RUN pip install .
