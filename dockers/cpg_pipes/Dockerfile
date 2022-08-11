FROM australia-southeast1-docker.pkg.dev/cpg-common/images/driver:4225655198684124178417573ffc6ed714cf0c04-hail-d4bd6a2fd606634057b0277fce06cedce5b764a6

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    gcc \
    build-essential \
    curl \
    && \
    add-apt-repository ppa:deadsnakes/ppa -y && \
    apt-get install -y --no-install-recommends python3.10 python3.10-distutils python3-apt python3.10-dev  && \
    curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10 && \
    pip3.10 install --no-cache-dir \
    cpg-pipes=0.3.8 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
