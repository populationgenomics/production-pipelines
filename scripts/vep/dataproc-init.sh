#!/bin/bash

##############
## Modified copy of https://github.com/hail-is/hail/blob/main/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh:
# * allow using a more recent VEP (e.g. 105), parametrised with $VEP_VERSION
# * using the CPG artifact registry VEP image, based on biocontainers image
# * changed the bucket to gs://cpg-reference/vep
##############

set -ex

export VEP_VERSION=__VEP_VERSION__
if [[ -z "$VEP_VERSION" ]]; then
    echo "Must provide VEP_VERSION in environment" 1>&2
    exit 1
fi

export IMAGE=australia-southeast1-docker.pkg.dev/cpg-common/images/vep:${VEP_VERSION}
export CONFIG_PATH=gs://cpg-reference/vep/${VEP_VERSION}/dataproc/config.json
export CACHE_PATH=gs://cpg-reference/vep/${VEP_VERSION}/cache.tar
export LOFTEE_PATH=gs://cpg-reference/vep/loftee.tar

mkdir -p /vep_data/loftee/

# Install docker
apt-get update
apt-get -y install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg2 \
    software-properties-common \
    tabix
curl -fsSL https://download.docker.com/linux/debian/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
apt-get update
apt-get install -y --allow-unauthenticated docker-ce

# Get VEP cache and LOFTEE data
gsutil cp ${CONFIG_PATH} /vep_data/vep-gcloud.json
# Will write /vep_data/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw
gsutil cat ${LOFTEE_PATH} | tar -xf - -C /vep_data/loftee/ &
# Will write /vep_data/vep/homo_sapiens/${VEP_VERSION}_GRCh38:
gsutil cat ${CACHE_PATH} | tar -xf - -C /vep_data/
# Copy the fasta to the top level so the path doesn't depend on VEP version
gunzip -c /vep_data/vep/homo_sapiens/*/Homo_sapiens.GRCh38*.fa.gz \
> /vep_data/Homo_sapiens.GRCh38.dna.fa

gcloud -q auth configure-docker australia-southeast1-docker.pkg.dev
docker pull ${IMAGE} &
wait

cat >/vep.c <<EOF
#include <unistd.h>
#include <stdio.h>
int
main(int argc, char *const argv[]) {
  if (setuid(geteuid()))
    perror( "setuid" );
  execv("/vep.sh", argv);
  return 0;
}
EOF
gcc -Wall -Werror -O2 /vep.c -o /vep
chmod u+s /vep

cat >/vep.sh <<EOF
#!/bin/bash
docker run -i -v /vep_data:/vep_data ${IMAGE} vep "\$@"
EOF
chmod +x /vep.sh

set +ex
