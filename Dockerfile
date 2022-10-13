FROM australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_utils:latest

COPY requirements.txt .
RUN pip install -r requirements.txt
