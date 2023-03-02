#!/usr/bin/env bash

# Sync the Broad reference resources into the corresponding CPG bucket.

gsutil rsync -r \
  gs://gcp-public-data--broad-references/hg38/v0/chrM/ \
  gs://cpg-common-main/references/hg38/v0/chrM
