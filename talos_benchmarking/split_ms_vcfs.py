"""
Oki, now we should have some multisample VCFs, split those out using bcftools +split

syntax template:
  - bcftools +split -Oz -W=tbi -o <folder> multisample.vcf.bgz
  - maybe shed a few reference sites with "bcftools view -c1"... but maybe don't

locally I needed to set a higher concurrent open-file limit
  - ulimit -n 10240
"""


import hail as hl

from cpg_utils import to_path, hail_batch


# two separate buckets, so use the seqr parent project
input_vcfs_1 = 'gs://cpg-ghfm-kidgen-test/talos_benchmarking/ghfm-kidgen_ms_vcfs'
input_vcfs_2 = 'gs://cpg-acute-care-test/talos_benchmarking/acute-care_ms_vcfs'

# write combined VCFs for both projects
output_vcfs = 'gs://cpg-seqr-test/talos_benchmarking/solo_vcfs'

batch_instance = hail_batch.get_batch('Generate single-sample VCFs from 2 projects')

image = 'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/talos:PR_552'

for proj, ms_vcf_source in [
    ('ghfm', input_vcfs_1),
    ('acute-care', input_vcfs_2),
]:

    all_files = to_path(ms_vcf_source).glob('*')

    for vcf_num, vcf_path in enumerate(all_files):
        this_job = batch_instance.new_job(f'Extracting chunk #{vcf_num} for proj {proj}')
        this_job.storage('200GiB').memory('32GiB').cpu(4)
        this_job.image(image)

        # localise the MS VCF
        ms_local = batch_instance.read_input(str(vcf_path))

        this_job.command(f"""
        mkdir $BATCH_TMPDIR/ss_outputs
        bcftools +split -Oz -W=tbi -o $BATCH_TMPDIR/ss_outputs {ms_local}
        
        gcloud storage cp "$BATCH_TMPDIR/ss_outputs/*" {output_vcfs}
        """)

batch_instance.run(wait=False)
