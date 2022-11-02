FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest

COPY README.md .
COPY setup.py .
COPY cpg_workflows cpg_workflows
RUN pip install .
