"""
Sync the Broad reference resources into the corresponding CPG bucket.
"""

gsutil rsync -r \
  gs://gatk-sv-resources-public/hg38/v0/sv-resources \
  gs://cpg-common-main/references/hg38/v0/sv-resources

gsutil cp \
  gs://broad-dsde-methods-eph/ped_1kgp_all.ped \
  gs://cpg-common-main/references/hg38/v0/sv-resources/ref-panel/
