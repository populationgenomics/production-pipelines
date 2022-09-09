""" 
Sync the Broad reference resources into the corresponding CPG bucket.
"""

gsutil rsync -r \
  gs://gatk-sv-resources-public/hg38/v0/sv-resources \
  gs://cpg-reference/hg38/v0/sv-resources
