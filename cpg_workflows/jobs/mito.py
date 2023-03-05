"""
Create Hail Batch jobs to run STRipy
"""
import hailtop.batch as hb
from hailtop.batch.job import Job
from textwrap import dedent

from cpg_utils import Path, to_path
from cpg_utils.hail_batch import image_path, fasta_res_group
from cpg_utils.hail_batch import command
from cpg_workflows.resources import STANDARD
from cpg_workflows.filetypes import CramPath
from cpg_workflows.utils import can_reuse
from cpg_workflows.targets import Sample

from cpg_workflows.jobs import picard


def subset_cram_to_chrM(
    b,
    cram_path: CramPath,
    output_bam_path: Path,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> tuple[Job | None, hb.ResourceGroup]:
    """
    Extract read pairs with at least one read mapping to chrM

    This uses the same extraction command as the broad pipeline
    as I am not clear the specifics *may* be important to how
    NUMT mapping reads are handled
    """

    # If we can resuse existing output, return the existing cram
    if can_reuse(
        [
            output_bam_path,
        ],
        overwrite,
    ):
        return None, b.read_input_group(
            **{
                'bam': str(output_bam_path),
                'bam.bai': str(output_bam_path) + '.bai',
            }
        )

    job_attrs = job_attrs or {}
    j = b.new_job('subset_cram_to_chrM', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    reference = fasta_res_group(b)
    j.declare_resource_group(
        output_bam={
            'bam': '{root}.bam',
            'bam.bai': '{root}.bai',
        }
    )

    cmd = f"""
        CRAM=$BATCH_TMPDIR/{cram_path.path.name}
        CRAI=$BATCH_TMPDIR/{cram_path.index_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path.path)} $CRAM
        retry_gs_cp {str(cram_path.index_path)} $CRAI

        gatk PrintReads \
            -R {reference.base} \
            -L chrM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -I $CRAM \
            --read-index $CRAI \
            -O {j.output_bam.bam}

        dirname {j.output_bam.bam}
        ls -l $(dirname {j.output_bam.bam})
    """

    j.command(command(cmd, define_retry_function=True))
    b.write_output(j.output_bam, str(output_bam_path.with_suffix('')))

    return j, j.output_bam


def mito_realign(
    b,
    sample_id: str,
    input_bam: hb.ResourceGroup,
    output_cram_path: Path,
    mito_ref: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job]:
    """
    Re-align reads to mito genome
    """
    if can_reuse(
        [
            output_cram_path,
        ],
        overwrite,
    ):
        return []

    job_attrs = job_attrs or {}
    j = b.new_job('mito_realign', job_attrs)
    j.image(image_path('bwa'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    nthreads = res.get_nthreads()

    cmd = f"""\
        bazam -Xmx16g \
            -n{min(nthreads, 6)} -bam {input_bam.bam} | \
        bwa \
            mem -K 100000000 -p -v 3 -t 2 -Y {mito_ref.base} \
            -R '@RG\\tID:{sample_id}\\tSM:{sample_id}' \
            - | \
        samtools view -bSu -T {mito_ref.base} - | \
        samtools sort -o {j.raw_cram}
        """
    j.command(command(cmd, define_retry_function=True))

    mkdup_j = picard.markdup(
        b=b,
        sorted_bam=j.raw_cram,
        output_path=output_cram_path,
        out_markdup_metrics_path=output_cram_path.with_suffix(
            '.markduplicates-metrics'
        ),
        fasta_reference=mito_ref
    )

    return [j, mkdup_j]


def mito_mutect2(
    b,
    cram_path: Path,
    reference: hb.ResourceGroup,
    region: str,
    max_reads_per_alignment_start: int = 75,
    job_attrs: dict | None = None,
) -> Job:
    """
    Call mutect2 on a mito genome

    """
    job_attrs = job_attrs or {}
    j = b.new_job('mito_mutect2', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)
    java_mem_mb = res.get_java_mem_mb()

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi', 'vcf.gz.stats': '{root}.vcf.gz.stats'}
    )

    cmd = f"""
        CRAM=$BATCH_TMPDIR/{cram_path.name}

        # Retrying copying to avoid google bandwidth limits
        retry_gs_cp {str(cram_path)} $CRAM
        retry_gs_cp {str(cram_path)}.crai $CRAM.crai

        # We need to create these files regardless, even if they stay empty
        # touch bamout.bam

        gatk --java-options "-Xmx{java_mem_mb}m" Mutect2 \
            -R {reference.base} \
            -I $CRAM \
            --read-filter MateOnSameContigOrNoMappedMateReadFilter \
            --read-filter MateUnmappedAndUnmappedReadFilter \
            -O {j.output_vcf['vcf.gz']} \
            --annotation StrandBiasBySample \
            --mitochondria-mode \
            --max-reads-per-alignment-start {max_reads_per_alignment_start} \
            --max-mnp-distance 0 \
            -L {region}
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def liftover_and_combine_vcfs(
    b,
    vcf: hb.ResourceGroup,
    shifted_vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    shift_back_chain: hb.ResourceFile,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Lifts over shifted vcf of control region and combines it with the rest of the chrM calls.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('liftover_and_combine_vcfs', job_attrs)
    j.image(image_path('picard'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        lifted_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        picard LiftoverVcf \
            I={shifted_vcf['vcf.gz']} \
            O={j.lifted_vcf['vcf.gz']} \
            R={reference.base} \
            CHAIN={shift_back_chain} \
            REJECT={j.rejected_vcf}.vcf.gz

        picard MergeVcfs \
            I={vcf['vcf.gz']} \
            I={j.lifted_vcf['vcf.gz']} \
            O={j.output_vcf['vcf.gz']}
    """

    j.command(command(cmd, define_retry_function=True))

    return j


def merge_mutect_stats(
    b,
    non_shifted_stats: hb.ResourceGroup,
    shifted_stats: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    Merge stats files from two mutect runs
    """
    job_attrs = job_attrs or {}
    j = b.new_job('merge_stats', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
            gatk MergeMutectStats \
                --stats {non_shifted_stats['vcf.gz.stats']} --stats {shifted_stats['vcf.gz.stats']} -O {j.combined_stats}

    """

    j.command(command(cmd, define_retry_function=True))

    return j


def filter_variants(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    merged_mutect_stats: hb.ResourceFile,
    max_alt_allele_count: int,
    vaf_filter_threshold: int,
    run_contamination: bool,
    contamination_estimate: hb.ResourceFile | None = None,
    f_score_beta: int | None = None,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    Mutect2 Filtering for calling Snps and Indels
    """
    job_attrs = job_attrs or {}
    j = b.new_job('filter_variants', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    blacklisted_sites = b.read_input_group(
        bed='gs://cpg-common-main/references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed',
        idx='gs://cpg-common-main/references/hg38/v0/chrM/blacklist_sites.hg38.chrM.bed.idx',
    )

    j.declare_resource_group(
        filtered_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    if run_contamination:
        contamination_estimate_string = f'--contamination-estimate  $(cat {contamination_estimate})'
    else:
        contamination_estimate_string = ''

    if f_score_beta is not None:
        f_score_beta_string = f'--f-score-beta  {f_score_beta}'
    else:
        f_score_beta_string = ''

    cmd = f"""
      gatk --java-options "-Xmx2500m" FilterMutectCalls -V {vcf['vcf.gz']} \
        -R {reference.base} \
        -O {j.filtered_vcf['vcf.gz']} \
        --stats {merged_mutect_stats} \
        --max-alt-allele-count {max_alt_allele_count} \
        --mitochondria-mode \
        --min-allele-fraction  {vaf_filter_threshold} \
        {contamination_estimate_string} {f_score_beta_string}

      gatk VariantFiltration -V {j.filtered_vcf['vcf.gz']} \
        -O {j.output_vcf['vcf.gz']} \
        --apply-allele-specific-filters \
        --mask {blacklisted_sites.bed} \
        --mask-name "blacklisted_site"

    """

    j.command(command(cmd, define_retry_function=True))

    return j


def split_multi_allelics_and_remove_non_pass_sites(
    b,
    vcf: hb.ResourceGroup,
    reference: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> Job:
    """
    split_multi_allelics_and_remove_non_pass_sites

    cmd taken from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L600
    """
    job_attrs = job_attrs or {}
    j = b.new_job('split_multi_allelics', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=4)
    res.set_to_job(j)

    j.declare_resource_group(
        split_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )
    j.declare_resource_group(
        output_vcf={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'}
    )

    cmd = f"""
        gatk LeftAlignAndTrimVariants \
            -R {reference.base} \
            -V {vcf['vcf.gz']} \
            -O {j.split_vcf['vcf.gz']} \
            --split-multi-allelics \
            --dont-trim-alleles \
            --keep-original-ac

        gatk SelectVariants \
            -V {vcf['vcf.gz']} \
            -O {j.output_vcf['vcf.gz']} \
            --exclude-filtered

        """

    j.command(command(cmd, define_retry_function=True))

    return j


def get_contamination(
    b,
    vcf: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    Uses new Haplochecker to estimate levels of contamination in mitochondria

    cmd taken from:
      https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L239

    CPG version returns results within files rather than strings supported in WDL.
    """
    job_attrs = job_attrs or {}
    j = b.new_job('get_contamination', job_attrs)
    # j.image(image_path('haplochecker'))
    j.image('us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124')

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""
        PARENT_DIR="$(dirname "{vcf['vcf.gz']}")"
        java -jar /haplocheckCLI/haplocheckCLI.jar "$PARENT_DIR"

        sed 's/\"//g' output > output-noquotes

        grep "SampleID" output-noquotes > headers
        FORMAT_ERROR="Bad contamination file format"
        if [ `awk '{{print $2}}' headers` != "Contamination" ]; then
            echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $6}}' headers` != "HgMajor" ]; then
            echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $8}}' headers` != "HgMinor" ]; then
            echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $14}}' headers` != "MeanHetLevelMajor" ]; then
            echo $FORMAT_ERROR; exit 1
        fi
        if [ `awk '{{print $15}}' headers` != "MeanHetLevelMinor" ]; then
            echo $FORMAT_ERROR; exit 1
        fi

        grep -v "SampleID" output-noquotes > output-data
        awk -F "\t" '{{print $2}}' output-data > {j.has_contamination}
        awk -F "\t" '{{print $6}}' output-data > {j.major_hg}
        awk -F "\t" '{{print $8}}' output-data > {j.minor_hg}
        awk -F "\t" '{{print $14}}' output-data > {j.major_level}
        awk -F "\t" '{{print $15}}' output-data > {j.minor_level}


        # Sanity check
        echo "has_contamination:"
        cat {j.has_contamination}

        echo "major_hg:"
        cat {j.major_hg}

        echo "minor_hg:"
        cat {j.minor_hg}

        echo "major_level:"
        cat {j.major_level}

        echo "minor_level:"
        cat {j.minor_level}

        # Instead of writhing to separate files lets munge it into json
        cat > {j.haplocheck_json} << EOF
        {{
            "has_contamination": "$(awk -F "\t" '{{print $2}}' output-data)",
            "major_hg": $(awk -F "\t" '{{print $6}}' output-data),
            "minor_hg": $(awk -F "\t" '{{print $8}}' output-data),
            "major_level": $(awk -F "\t" '{{print $14}}' output-data),
            "minor_level": $(awk -F "\t" '{{print $15}}' output-data)
        }}
        EOF
        """

    j.command(command(cmd, define_retry_function=True))

    return j

def get_contamination_dummy(
    b,
    vcf: hb.ResourceGroup,
    job_attrs: dict | None = None,
) -> Job:
    """
    dumy job for testing
    """
    job_attrs = job_attrs or {}
    j = b.new_job('get_contamination_dummy', job_attrs)
    j.image(image_path('gatk'))

    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    cmd = f"""

        # Instead of writhing to separate files lets munge it into json
        cat > {j.haplocheck_json} << EOF
        {{
            "has_contamination": "YES",
            "major_hg": "foo",
            "minor_hg": "bar",
            "major_level": 0.05,
            "minor_level": 0.03
        }}
        EOF
        """

    j.command(command(cmd, define_retry_function=True))

    return j

def _get_max_contamination(
    haplocheker_json: hb.ResourceFile,
    # verify_bam_id_result: hb.ResourceGroup,
) -> float:
    """
    Parse the output from haplocheker to define the estimated contamination level
    to use for filtering.

    Implements the logic from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L523-L525
    """
    import json

    with open(str(haplocheker_json)) as f:
        haplocheker = json.load(f)

        if  haplocheker['has_contamination'] == "YES":
            if haplocheker['contamination_major'] == 0.0:
                max_contamination = haplocheker['contamination_minor']
            else:
                max_contamination = 1.0 - haplocheker['contamination_major']
        else:
            max_contamination = 0.0

        # TODO: look at verify_bam_id and integrate that as well
        # if verify_bam_id > max_contamination:
        #     max_contamination = verify_bam_id

        return max_contamination


def get_max_contamination(
    b,
    haplocheker_json: hb.ResourceFile,
    # verify_bam_id_result: hb.ResourceFile,
    job_attrs: dict | None = None,
) -> tuple[Job, hb.ResourceFile]:
    """
    Parse the output from haplocheker to define the estimated contamination level
    to use for filtering.

    It should also incorporate the estimate from verify_bam_id if it is available... but it does not yet (TODO).

    Implements the logic from:
    https://github.com/broadinstitute/gatk/blob/227bbca4d6cf41dbc61f605ff4a4b49fc3dbc337/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L523-L525
    """

    j = b.new_python_job('get_max_contamination_j', job_attrs)
    j.image(image_path('cpg_workflows'))
    res = STANDARD.request_resources(ncpu=2)
    res.set_to_job(j)

    contamination_estimate = j.call(get_max_contamination, haplocheker_json)

    return j, contamination_estimate.as_json


def genotype_mito(
    b,
    cram_path: Path,
    shifted_cram_path: Path,
    output_vcf_path: Path,
    # haplochecker_json_path: Path,
    mito_reff: hb.ResourceGroup,
    shifted_mito_reff: hb.ResourceGroup,
    job_attrs: dict | None = None,
    overwrite: bool = False,
) -> list[Job | None]:
    """
    Genotype mitochondrial genome using mutect2

    outputs:
     - output_vcf_path: Final vcf for individual.
     - haplochecker_json_path: Contamination levels determined by haplochecker tool.


    Re implementation of the functions from:
      https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L89

    This set of jobs:
        - Calls variants for all regions of chrM except the control region from the cram mapped to the normal reference.
        - Calls variants for in only the control region from the cram mapped to the shifted reference.
        - Uses LiftOver to lift the control region variants back to normal chrM coordinate space.
        - Merges the two vcfs into a single vcf.
        - Merges the Mutect stats output files from both call sets.
        - Performs an initial round of variant filtering.
        - Uses haplotype based analysis to estimate level of MT contamination.
        - Performs an second round of variant filtering based on the estimated contamination.
    """

    # Call variants on WT genome
    call_j = mito_mutect2(
        b=b,
        cram_path=cram_path,
        reference=mito_reff,
        region='chrM:576-16024',  # Exclude the control region.
        job_attrs=job_attrs,
    )

    # Call variants in the control region using a shifted reference
    shifted_call_j = mito_mutect2(
        b=b,
        cram_path=shifted_cram_path,
        reference=shifted_mito_reff,
        region='chrM:8025-9144',  # Only call inside the control region.
        job_attrs=job_attrs,
    )
    jobs = [call_j, shifted_call_j]

    # Merge the wt and shifted VCFs
    merge_j = liftover_and_combine_vcfs(
        b=b,
        vcf=call_j.output_vcf,
        shifted_vcf=shifted_call_j.output_vcf,
        reference=mito_reff,  # Does this need to be full genome ref?
        shift_back_chain=shifted_mito_reff.shift_back_chain,
        job_attrs=job_attrs,
    )
    jobs.append(merge_j)

    # Merge the mutect stats output files (needed for filtering)
    merge_stats_J = merge_mutect_stats(
        b=b,
        non_shifted_stats=call_j.output_vcf,
        shifted_stats=shifted_call_j.output_vcf,
    )
    jobs.append(merge_stats_J)

    # Initial round of filtering to exclude blacklist and high alt alleles
    initial_filter_j = filter_variants(
        b=b,
        vcf=merge_j.output_vcf,
        reference=mito_reff,
        merged_mutect_stats=merge_stats_J.combined_stats,
        # alt_allele and vaf config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
        max_alt_allele_count=4,
        vaf_filter_threshold=0,
        run_contamination=False,
        job_attrs=job_attrs
    )
    jobs.append(initial_filter_j)

    # SplitMultiAllelicsAndRemoveNonPassSites
    split_multiallelics_j = split_multi_allelics_and_remove_non_pass_sites(
        b=b,
        vcf=initial_filter_j.output_vcf,
        reference=mito_reff,
    )
    jobs.append(split_multiallelics_j)

    # Use mito reads to identify level of contamination
    # get_contamination_j = get_contamination(
    #     b=b,
    #     vcf=split_multiallelics_j.output_vcf,
    # )
    # jobs.append(get_contamination_j)
    get_contamination_j = get_contamination_dummy(  ######### LOOK HERE!!!!!!!!!!!
        b=b,
        vcf=split_multiallelics_j.output_vcf,
    )
    jobs.append(get_contamination_j)

    # Use haplochecker results to determine final contamination filtering threshold
    get_max_contamination_j,  contamination_estimate = get_max_contamination(
        b=b,
        haplocheker_json=get_contamination_j.haplocheck_json
    )

    # b.new_python_job('get_max_contamination_j', job_attrs)
    # contamination_estimate = get_max_contamination_j.call(get_max_contamination, get_contamination_j.haplocheck_json)
    # jobs.append(get_max_contamination_j)

    # Filter round 2 - remove contamination
    second_filter_j = filter_variants(
        b=b,
        vcf=initial_filter_j.output_vcf,
        reference=mito_reff,
        merged_mutect_stats=merge_stats_J.combined_stats,
        # alt_allele and vaf config from https://github.com/broadinstitute/gatk/blob/master/scripts/mitochondria_m2_wdl/AlignAndCall.wdl#L167
        max_alt_allele_count=4,
        vaf_filter_threshold=0,
        run_contamination=True,
        contamination_estimate=contamination_estimate,
        job_attrs=job_attrs
    )
    jobs.append(second_filter_j)

    return jobs
