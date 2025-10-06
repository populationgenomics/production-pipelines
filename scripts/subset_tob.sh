#!/bin/bash

gcloud storage ls gs://cpg-tob-wgs-main/ica/dragen_3_7_8/output/recal_gvcf/*.gz > gvcf_files

while read line
do
    new_name=$(echo $line | rev | cut -d / -f 1 | rev | sed 's/hard/chrM.hard/')
    bcftools view -r chrM -o "$new_name" -Oz "$line"
    rm *.tbi
    bcftools index "$new_name"
done < gvcf_files

gcloud storage cp *chrM.hard* gs://cpg-tob-wgs-main/ica/dragen_3_7_8/output/recal_gvcf/chrm/
