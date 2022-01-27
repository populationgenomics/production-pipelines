#!/bin/bash

##############
## Modified copy of https://github.com/hail-is/hail/blob/main/hail/python/hailtop/hailctl/dataproc/resources/vep-GRCh38.sh:
# * using a more recent VEP (105 vs 95)
# * using cpg artifact registry VEP image, based on biocontainers image
# * changed hardcoded hail bucket to cpg-reference/vep
##############

export PROJECT="$(gcloud config get-value project)"
export VEP_CONFIG_PATH=/vep_data/vep-gcloud.json
export VEP_BUCKET=cpg-reference/vep
export ASSEMBLY=GRCh38
export VEP_DOCKER_IMAGE=australia-southeast1-docker.pkg.dev/cpg-common/images/vep:105

mkdir -p /vep_data/loftee_data
mkdir -p /vep_data/homo_sapiens

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
gsutil -u $PROJECT cp gs://${VEP_BUCKET}/vep95-GRCh38-loftee-gcloud.json /vep_data/vep105-GRCh38-gcloud.json
ln -s /vep_data/vep105-GRCh38-gcloud.json $VEP_CONFIG_PATH

gsutil -u $PROJECT cat gs://${VEP_BUCKET}/loftee-beta.tar | tar -xf - -C /vep_data/ &
gsutil -u $PROJECT cat gs://${VEP_BUCKET}/homo_sapiens_vep_105_GRCh38.tar.gz | tar -xf - -C /vep_data/homo_sapiens
ls /vep_data/
ls /opt/vep/.vep/homo_sapiens
ls /opt/vep/.vep/homo_sapiens/105_GRCh38

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
docker run -i -v /vep_data/:/opt/vep/.vep/:ro ${VEP_DOCKER_IMAGE} \
  /usr/local/bin/vep "\$@"
EOF
chmod +x /vep.sh
