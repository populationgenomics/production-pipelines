#!/usr/bin/env bash

gcloud storage cp --recursive \
    gs://gcp-public-data--broad-references/hg38/v0/chrM \
    gs://cpg-common-main/references/hg38/v0/chrM
