"""
Tests to reproduce issues reported in https://github.com/populationgenomics/seqr-private/issues/1
"""
from cpg_pipes import Namespace, to_path, images
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs import split_intervals, align
from cpg_pipes.jobs.align import Aligner, MarkDupTool
from cpg_pipes.jobs.haplotype_caller import produce_gvcf, merge_gvcfs_job
from cpg_pipes.jobs.somalier import pedigree
from cpg_pipes.pipeline import create_pipeline
from cpg_pipes.targets import Dataset, Sample
from cpg_pipes.types import CramPath, AlignmentInput, FastqPair

SID = 'CPG54353'
# SID = 'CPG54650'
# SID = 'CPG54130'

BUCKET = to_path('gs://cpg-fewgenomes-test/tmp/acute-care-main')
UPLOAD_BUCKET = to_path('gs://cpg-fewgenomes-test/tmp/acute-care-main-upload')


def main():
    pipe = create_pipeline(
        analysis_dataset='fewgenomes',
        name=f'test-missing-variant-{SID}',
        version='v0',
        namespace=Namespace.TEST,
        description=f'Test missing variant in {SID}',
    )

    # cram_path = CramPath(BUCKET / f'cram/bwamem/{SID}.cram')
    cram_path = CramPath(BUCKET / f'cram/nagim/{SID}.cram')

    # intervals_j, intervals = split_intervals.get_intervals(
    #     b=pipe.b,
    #     refs=pipe.refs,
    #     scatter_count=50,
    # )
    # 
    # out_gvcf = to_path(bucket / 'gvcf-v4230-dragen-false' / 'out.g.vcf.gz')
    # jobs4230 = produce_gvcf(
    #     pipe.b,
    #     sample_name=f'{sid}-v4230-dragen-false',
    #     refs=pipe.refs,
    #     cram_path=cram_path,
    #     intervals=intervals,
    #     scatter_count=50,
    #     tmp_bucket=bucket / 'gvcf-v4230-dragen-false' / 'tmp',
    #     output_path=out_gvcf,
    #     dragen_mode=False,
    # )
    # 
    # out_gvcf = to_path(bucket / 'gvcf-v4230-dragen-true' / 'out.g.vcf.gz')
    # jobs4230_dragen = produce_gvcf(
    #     pipe.b,
    #     sample_name=f'{sid}-v4230-dragen-true',
    #     refs=pipe.refs,
    #     cram_path=cram_path,
    #     intervals=intervals,
    #     scatter_count=50,
    #     tmp_bucket=bucket / 'gvcf-v4230-dragen-true' / 'tmp',
    #     output_path=out_gvcf,
    #     dragen_mode=True,
    # )
    
    # gatk_reproduce_exactly_job(
    #     pipe,
    #     cram_path,
    #     sample_name=f'{SID}-bwamem-v4210-dragen-false-reproduce-exactly',
    #     out_gvcf_path=BUCKET / f'{SID}-v4210-dragen-false-reproduce-exactly' / 'out.g.vcf.gz'
    # )
    
    # re-alignment?
    # TODO: start a docker and only re-align the region?
    alignment_input = [
        FastqPair(
            UPLOAD_BUCKET / 'cpg_acute_positives_20211003_213917/200212_A00692_0063_ML201946_20W000180-FAM000358_MAN-20200212_NEXTERAFLEXWGS_L001_R1.fastq.gz',
            UPLOAD_BUCKET / 'cpg_acute_positives_20211003_213917/200212_A00692_0063_ML201946_20W000180-FAM000358_MAN-20200212_NEXTERAFLEXWGS_L001_R2.fastq.gz',
        ),
        FastqPair(
            UPLOAD_BUCKET / 'cpg_acute_positives_20211003_213917/200212_A00692_0063_ML201946_20W000180-FAM000358_MAN-20200212_NEXTERAFLEXWGS_L002_R1.fastq.gz',
            UPLOAD_BUCKET / 'cpg_acute_positives_20211003_213917/200212_A00692_0063_ML201946_20W000180-FAM000358_MAN-20200212_NEXTERAFLEXWGS_L002_R2.fastq.gz',
        ),
    ]
    sid = 'CPG54353'
    # align.align(
    #     pipe.b,
    #     alignment_input=alignment_input,
    #     sample_name=sid,
    #     refs=pipe.refs,
    #     output_path=BUCKET / 'realigned' / f'{sid}-bwa-biobambam.cram',
    # )
    # align.align(
    #     pipe.b,
    #     alignment_input=alignment_input,
    #     sample_name=sid,
    #     refs=pipe.refs,
    #     aligner=Aligner.DRAGMAP,
    #     markdup_tool=MarkDupTool.BIOBAMBAM,
    #     output_path=BUCKET / 'realigned' / f'{sid}-dragmap-picard.cram',
    # )
    
    _index_cram(pipe.b, BUCKET / 'cram/bwamem/CPG54353.cram')

    # _gatk_reproduce_exactly_job_hc(
    #     pipe=pipe,
    #     cram_path=cram_path,
    #     sample_name=SID,
    #     interval=pipe.b.read_input('gs://cpg-fewgenomes-test-tmp/test-missing-variant-CPG54353/v0/hail/26cefd/Make_50_intervals-m68kK/intervals/0045-scattered.interval_list'),
    #     out_gvcf_path=BUCKET / f'{SID}-hc-0045-nagimcram.g.vcf.gz',
    #     # localise_cram=False,  # not localising the CRAM reproduces the issue
    #     localise_cram=False,  # not localising the CRAM reproduces the issue
    # )
    
    # run_somalier(pipe)

    pipe.run(keep_scratch=True)


def _index_cram(b, cram):
    j = b.new_job('Re-index')
    j.storage('50G')
    j.cpu(16)
    j.image(images.SAMTOOLS_PICARD_IMAGE)
    cram = b.read_input(str(cram))
    j.command(f"""
    samtools index {cram}
    cp {cram}.crai {j.out_crai}
    """)
    b.write_output(j.out_crai, str(cram) + '.crai')


def run_somalier(pipe):
    ds = Dataset.create(
        'test-missing-variant', 
        namespace=Namespace.TEST,
    )
    gvcf_by_sid = {
        # 'CPG54353-cram':          BUCKET / 'cram/CPG54353.cram',
        # 'CPG54353-cram-nagim':    BUCKET / 'cram/nagim/CPG54353.cram',
        # 'CPG54353-cram-bwamem':   BUCKET / 'cram/bwamem/CPG54353.cram',
        'CPG54353-gvcf':          BUCKET / 'gvcf/CPG54353.somalier',
        'CPG54353-gvcf-nagim':    BUCKET / 'gvcf/nagim/CPG54353.g.vcf.gz',
        'CPG54353-gvcf-gatk4210': BUCKET / 'gvcf/gatk-4210/CPG54353.g.vcf.gz',
        'CPG54353-reprocessed':   BUCKET / 'CPG54353-v4210-dragen-false-reproduce-exactly/out.g.vcf.gz',
        'CPG54650-gvcf':          BUCKET / 'gvcf/CPG54650.g.vcf.gz',
        'CPG54650-gvcf-nagim':    BUCKET / 'gvcf/nagim/CPG54650.g.vcf.gz',
        'CPG54650-gvcf-gatk4210': BUCKET / 'gvcf/gatk-4210/CPG54650.g.vcf.gz',
        'CPG54650-reprocessed':   BUCKET / 'CPG54650-v4210-dragen-false-reproduce-exactly/out.g.vcf.gz',
    }
    for k, v in gvcf_by_sid.items():
        ds.add_sample(k)

    pedigree(
        b=pipe.b,
        dataset=ds,
        input_path_by_sid=gvcf_by_sid,
        refs=pipe.refs,
        overwrite=True,
        out_samples_path=BUCKET / 'somalier-samples.tsv',
        out_pairs_path=BUCKET / 'somalier-pairs.tsv',
        out_html_path=BUCKET / 'somalier.html',
        tmp_bucket=BUCKET / 'tmp',
    )


def gatk_reproduce_exactly_job(
    pipe,
    cram_path,
    sample_name,
    out_gvcf_path,
):
    intervals_j = _add_split_intervals_job(
        b=pipe.b,
        interval_list='gs://cpg-reference/hg38/v0/hg38.even.handcurated.20k.intervals',
        scatter_count=50,
        ref_fasta=pipe.refs.ref_fasta,
    )
    
    gvcfs = []
    hc_jobs = []
    for idx in range(50):
        idx_raw_gvcf_path = str(out_gvcf_path).replace('.g.vcf.gz', f'-{idx}-raw.g.vcf.gz')
        j1 = _gatk_reproduce_exactly_job_hc(
            pipe=pipe,
            cram_path=cram_path,
            sample_name=sample_name,
            interval=intervals_j.intervals[f'interval_{idx}'],
            out_gvcf_path=idx_raw_gvcf_path,
        )
        hc_jobs.append(j1)
        gvcfs.append(j1.output_gvcf)
    merged_gvcf_path = str(out_gvcf_path).replace('.g.vcf.gz', f'-merged.g.vcf.gz')
    merge_j = merge_gvcfs_job(
        b=pipe.b,
        sample_name=sample_name,
        gvcfs=gvcfs,
        out_gvcf_path=to_path(merged_gvcf_path),
        overwrite=True,
    )
    merge_j.depends_on(*hc_jobs)
    postproc_j = _gatk_reproduce_exactly_job_reblock(
        pipe=pipe,
        raw_gvcf=merge_j.output_gvcf['g.vcf.gz'],
        sample_name=sample_name,
        out_gvcf_path=out_gvcf_path,
    )
    postproc_j.depends_on(merge_j)
    # gvcfs.append(j2.output_gvcf)
    # gvcfs.append(pipe.b.read_input_group(**{
    #     'g.vcf.gz': str(idx_gvcf_path),
    #     'g.vcf.gz.tbi': str(idx_gvcf_path) + '.tbi',
    # }))

    return merge_j


def _add_split_intervals_job(
    b,
    interval_list: str,
    scatter_count: int,
    ref_fasta: str,
):
    """
    Split genome into intervals to parallelise GnarlyGenotyper.
    Returns: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f'Make {scatter_count} intervals')
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:4.2.1.0')
    java_mem = 3
    j.memory('standard')  # ~ 4G/core ~ 4G
    j.storage('16G')
    j.declare_resource_group(
        intervals={
            f'interval_{idx}': f'{{root}}/'
                               f'{str(idx).zfill(4)}-scattered.interval_list'
            for idx in range(scatter_count)
        }
    )

    j.command(
        f"""set -e
    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number 
    # of intervals. But we have to expect exactly the {scatter_count} number of 
    # output files because our workflow is not dynamic.
    gatk --java-options -Xms{java_mem}g SplitIntervals \\
      -L {interval_list} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {ref_fasta} \\
      -mode INTERVAL_SUBDIVISION
      """
    )
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


def _gatk_reproduce_exactly_job_hc(
    pipe,
    cram_path,
    sample_name,
    interval,
    out_gvcf_path,
    localise_cram=False,
):
    j = pipe.b.new_job(f'{sample_name}: gatk 4210 dragen false: reproduce exactly (HaplotypeCaller)')
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:4.2.1.0')
    res = STANDARD.set_resources(j, storage_gb=45)
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )

    # CRAM=/io/batch/{sample_name}.cram
    # CRAI=/io/batch/{sample_name}.cram.crai
    # # Retrying copying to avoid google bandwidth limits
    # retry_gs_cp {str(cram_path.path)} $CRAM
    # retry_gs_cp {str(cram_path.index_path)} $CRAI
    # --output {j.output_gvcf['g.vcf.gz']} --intervals {interval} --input {cram_path} --read-index {cram_path}.crai --reference 

    if localise_cram:
        cram = pipe.b.read_input(str(cram_path.path))
        crai = pipe.b.read_input(str(cram_path.path) + '.crai')
    else:
        cram = str(cram_path)
        crai = str(cram_path) + '.crai'

    reference = pipe.refs.fasta_res_group(pipe.b)
    cmd = f"""\
    export GOOGLE_APPLICATION_CREDENTIALS=/gsa-key/key.json
    gcloud -q auth activate-service-account --key-file=$GOOGLE_APPLICATION_CREDENTIALS

    gatk --java-options "-Xms7g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
    HaplotypeCaller \\
    -R {reference.base} \\
    -I {cram} \\
    --read-index {crai} \\
    -L {interval} \\
    --disable-spanning-event-genotyping \\
    -O {j.output_gvcf['g.vcf.gz']} \\
    -G AS_StandardAnnotation \\
    -GQB 20 \\
    -ERC GVCF
    """
    j.command(wrap_command(cmd))
    pipe.b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j
   
   
def _gatk_reproduce_exactly_job_reblock(
    pipe,
    raw_gvcf,
    sample_name,
    out_gvcf_path,
):
    j = pipe.b.new_job(f'{sample_name}: gatk 4210 dragen false: reproduce exactly (ReblockGVCF)')
    j.image('australia-southeast1-docker.pkg.dev/cpg-common/images/gatk:4.2.1.0')
    res = STANDARD.set_resources(j, storage_gb=45)
    j.declare_resource_group(
        output_gvcf={
            'g.vcf.gz': '{root}-' + sample_name + '.g.vcf.gz',
            'g.vcf.gz.tbi': '{root}-' + sample_name + '.g.vcf.gz.tbi',
        }
    )
    reference = pipe.refs.fasta_res_group(pipe.b)
    noalt_regions = pipe.b.read_input(str(pipe.refs.noalt_regions))
    cmd = f"""
    GVCF={raw_gvcf}
    REBLOCKED=/io/batch/{sample_name}-reblocked.g.vcf.gz

    tabix -f -p vcf $GVCF

    gatk --java-options -Xms7g ReblockGVCF --reference {reference.base} -V $GVCF -do-qual-approx -O $REBLOCKED --create-output-variant-index true
    
    cp $REBLOCKED  {j.output_gvcf['g.vcf.gz']}

    tabix -p vcf {j.output_gvcf['g.vcf.gz']}
    """
    j.command(
        wrap_command(
            cmd, setup_gcp=True, monitor_space=True, define_retry_function=True
        )
    )
    pipe.b.write_output(j.output_gvcf, str(out_gvcf_path).replace('.g.vcf.gz', ''))
    return j


main()
