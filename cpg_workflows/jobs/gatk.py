import hailtop.batch as hb
from hailtop.batch.job import Job

from cpg_utils import Path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import command, fasta_res_group, image_path
from cpg_workflows.utils import can_reuse


def preprocess_intervals(b, intervals_path, job_attrs, output_path):
    if can_reuse(output_path):
        return []

    job_attrs = (job_attrs or {}) | {
        'tool': 'gatk PreprocessIntervals',
    }

    j = b.new_job('Preprocess intervals', job_attrs)
    j.image(image_path('gatk'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(b)

    if sequencing_type == 'exome':
        cmd = f"""
        gatk PreprocessIntervals --reference {reference.base} --intervals {intervals_path} \\
            --padding 250 --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY \\
            --output {j.preprocessed_intervals}
        """

    elif sequencing_type == 'genome':
        cmd = f"""
        gatk PreprocessIntervals --reference {reference.base} \\
            --interval-merging-rule OVERLAPPING_ONLY --padding 0 \\
            --output {j.preprocessed_intervals}
        """

    else:
        raise ValueError('workflow/sequencing_type must be one of: "exome", "genome"')

    j.command(command(cmd))
    b.write_output(j.preprocessed_intervals, output_path)
    return [j]


def collect_read_counts(b, sample, intervals, job_attrs, output_path):
    if can_reuse(output_path):
        return []

    job_attrs = (job_attrs or {}) | {
        'tool': 'gatk CollectReadCounts',
    }

    j = b.new_job(f'Collect gCNV counts for {sample.id}', job_attrs)
    j.image(image_path('gatk'))

    reference = fasta_res_group(b)
    cram_path = sample.make_cram_path()
    cmd = f"""
    gatk CollectReadCounts \\
      --reference {reference.base} --intervals {intervals} \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --input {str(cram_path.path)} --read-index {str(cram_path.index_path)} \\
      --output {j.gcnv_counts}
    """

    j.command(command(cmd, setup_gcp=True))
    b.write_output(j.gcnv_counts, output_path)
    return [j]


def determine_ploidy(b, cohort_name, ploidy_priors, inputs, job_attrs, output_dir):
    job_attrs = (job_attrs or {}) | {
        'tool': 'gatk DetermineGermlineContigPloidy',
    }

    j = b.new_job('Determine ploidy for {cohort_name}', job_attrs)
    j.image(image_path('gatk'))

    input_args = ' '.join([f'--input {f}' for f in inputs])
    cmd = f"""
    gatk DetermineGermlineContigPloidy \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --contig-ploidy-priors {ploidy_priors} \\
      {input_args} \\
      --output $BATCH_TMPDIR --output-prefix {cohort_name}
    """
    # TODO maybe add --intervals {intervals_path}
    #   tar czf ~{cohort_name}-contig-ploidy-model.tar.gz -C ~{output_dir_}/~{cohort_name}-model .
    #   tar czf ~{cohort_name}-contig-ploidy-calls.tar.gz -C ~{output_dir_}/~{cohort_name}-calls .

    j.command(command(cmd, setup_gcp=True))
    return [j]
