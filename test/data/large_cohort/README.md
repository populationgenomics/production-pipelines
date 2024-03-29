To prepare data:

```bash
# Prepare GVCFs
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
mkdir -p gvcf
for s in SAMPLE1 SAMPLE2
do
  bcftools view -R intervals.bed gs://cpg-thousand-genomes-main/gvcf/$s.g.vcf.gz -Oz -o gvcf/$s.g.vcf.gz
  tabix gvcf/$s.g.vcf.gz
done

# Prepare reference
python build_ref_tables.py
gsutil cp -r gs://cpg-common-main/references/gnomad/v0-toy-chr20-x-y reference/gnomad/v0

mkdir -p reference/hg38/v0
gsutil cp gs://cpg-common-main/references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  tmp-wgs_calling_regions.hg38.interval_list
picard BedToIntervalList -I intervals.bed \
  --SEQUENCE_DICTIONARY /Users/vlad/bio/hg38/Homo_sapiens_assembly38.dict \
  -O intervals.interval_list
picard IntervalListTools --ACTION INTERSECT \
  -I tmp-wgs_calling_regions.hg38.interval_list \
  -I intervals.interval_list \
  -O reference/hg38/v0/wgs_calling_regions.hg38.interval_list

GTF=gs://cpg-common-main/references/hg38/v0/sv-resources/resources/v1/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf
mkdir -p reference/hg38/v0/sv-resources/resources/v1
bedtools intersect -wa \
  -a <(gsutil cat $GTF) \
  -b intervals.bed > \
  reference/hg38/v0/sv-resources/resources/v1/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf
```
