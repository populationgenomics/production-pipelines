#!/usr/bin/env bash

set -ex

gcloud auth configure-docker australia-southeast1-docker.pkg.dev
skopeo copy docker://us.gcr.io/broad-dsde-methods/markw/gatk:2023-07-13-4.4.0.0-43-gd79823f9c-NIGHTLY-SNAPSHOT docker://australia-southeast1-docker.pkg.dev/cpg-common/images/sv/gatk:2023-07-13-4.4.0.0-43-gd79823f9c-NIGHTLY-SNAPSHOT