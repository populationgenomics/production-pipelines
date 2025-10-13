#!/usr/bin/env bash

# script for using gcloud to localise and delocalise a file (test of multipart upload performance/timing)
# unlike hail batch this isn't timed within the interface, so we require timestamping in the script

in_file=$1
out_file=$2

base_name=$(basename "${in_file}")

gcloud config set storage/parallel_composite_upload_enabled False

gcloud storage cp "${in_file}" "${BATCH_TMPDIR}/${base_name}"

time gcloud storage cp "${BATCH_TMPDIR}/${base_name}" "${out_file}"
