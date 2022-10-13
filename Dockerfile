FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest

COPY requirements.txt .
RUN python3 -m pip install -r requirements.txt
