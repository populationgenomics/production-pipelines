from os.path import basename

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
    j.image(image_path('gatk_gcnv'))

    sequencing_type = get_config()['workflow']['sequencing_type']
    reference = fasta_res_group(b)

    exclude_intervals = ' '.join([f'--exclude-intervals {i}' for i in get_config()['workflow'].get('exclude_intervals', [])])

    if sequencing_type == 'exome':
        cmd = f"""
        gatk PreprocessIntervals \\
          --reference {reference.base} --intervals {intervals_path} {exclude_intervals} \\
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
    b.write_output(j.preprocessed_intervals, str(output_path))
    return [j]


def collect_read_counts(b, sample, intervals, job_attrs, output_path):
    if can_reuse(output_path):
        return []

    job_attrs = (job_attrs or {}) | {
        'tool': 'gatk CollectReadCounts',
    }

    j = b.new_job(f'Collect gCNV counts for {sample.id}', job_attrs)
    j.image(image_path('gatk_gcnv'))

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
    b.write_output(j.gcnv_counts, str(output_path))
    return [j]


def determine_ploidy(b: hb.Batch, cohort_name: str, ploidy_priors: str, inputs: list[Path], job_attrs: dict, output_dir):
    """

    Args:
        b (hb.Batch):
        cohort_name (str):
        ploidy_priors (str): file from config
        inputs (list[Path]): CollectReadCounts results for each Sample
        job_attrs (dict):
        output_dir (Path):

    Returns:

    """
    job_attrs = (job_attrs or {}) | {
        'tool': 'gatk DetermineGermlineContigPloidy',
    }

    j = b.new_job(f'Determine ploidy for {cohort_name}', job_attrs)
    j.image(image_path('gatk_gcnv'))

    input_args = ''
    for f in inputs:
        sample_input = b.read_input(str(f))
        input_args += f' --input {sample_input}'

    # --contig-ploidy-priors argument must be a local file
    ploidy_prior = b.read_input(ploidy_priors)

    # define output locations in advance
    j.declare_resource_group(
        contig_model={
            '-calls.tar.gz': '{root}-calls.tar.gz',
            '-model.tar.gz': '{root}-model.tar.gz'
        }
    )

    cmd = f"""
    gatk DetermineGermlineContigPloidy \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --contig-ploidy-priors {ploidy_prior} \\
      {input_args} \\
      --output $BATCH_TMPDIR --output-prefix {j.contig_model}
    """

    # writes outputs to cloud location, subbing in {root} for output_dir
    b.write_output(j.contig_model, str(output_dir))

    # TODO maybe add --intervals {intervals_path}
    #   tar czf ~{cohort_name}-contig-ploidy-model.tar.gz -C ~{output_dir_}/~{cohort_name}-model .
    #   tar czf ~{cohort_name}-contig-ploidy-calls.tar.gz -C ~{output_dir_}/~{cohort_name}-calls .

    j.command(command(cmd, setup_gcp=True, define_retry_function=True))
    return [j]
