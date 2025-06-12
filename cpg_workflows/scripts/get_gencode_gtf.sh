#!/usr/bin/env bash

# use this get_gencode_gtf.sh script to download the latest GENCODE GTF file
# and place it in the appropriate place in the CPG bucket

set -ex

GENCODE_VERSION=${1:-"47"}

# download the GTF file
wget -O gencode.gtf.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz"

# don't unzip the file - unzipping inflates from ~50MB to ~1.8GB, and transfer is probably the bottleneck
# place the file in the place in the CPG bucket
gcloud storage cp gencode.gtf.gz "gs://cpg-common-main/references/hg38/v0/gencode_${GENCODE_VERSION}.gtf.gz"
