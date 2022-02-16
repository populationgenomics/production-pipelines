#!/bin/bash

##############
## Modified copy of https://github.com/hail-is/hail/blob/main/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh:
# * using a more recent VEP (105 vs. 95)
# * using cpg artifact registry VEP image, based on biocontainers image
# * changed hardcoded hail bucket to cpg-reference/vep
##############

export PROJECT="$(gcloud config get-value project)"
export VEP_CONFIG_PATH=/vep_data/vep-gcloud.json
export VEP_BUCKET=cpg-reference/vep
export ASSEMBLY=GRCh38
export VEP_DOCKER_IMAGE=quay.io/biocontainers/ensembl-vep:105.0--pl5262h4a94de4_0

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
gsutil -u $PROJECT cp gs://${VEP_BUCKET}/vep105-GRCh38-loftee-gcloud.json $VEP_CONFIG_PATH
# Will write /vep_data/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw
gsutil -u $PROJECT cat gs://${VEP_BUCKET}/loftee-beta.tar | tar -xf - -C /vep_data/loftee/ &
# Will write /vep_data/vep/homo_sapiens/105_GRCh38:
gsutil -u $PROJECT cat gs://${VEP_BUCKET}/homo_sapiens_vep_105_GRCh38.tar | tar -xf - -C /vep_data/

docker pull ${VEP_DOCKER_IMAGE} &
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
docker run -i \
-v /vep_data:/vep_data:ro \
-v $VEP_CONFIG_PATH:$VEP_CONFIG_PATH:ro \
${VEP_DOCKER_IMAGE} /usr/local/bin/vep "\$@"
EOF
chmod +x /vep.sh
