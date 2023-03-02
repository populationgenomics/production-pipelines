"""
Sync the Broad reference resources into the corresponding CPG bucket.
"""

gsutil rsync -r \
  gs://gatk-sv-resources-public/hg38/v0/chrM \
  gs://cpg-common-main/references/hg38/v0/chrM
