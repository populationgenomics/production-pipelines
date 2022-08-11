"""
Adding jobs for fingerprinting and pedigree checks. Mostly using Somalier.
"""
import logging
from os.path import basename

import pandas as pd
from cpg_utils.hail_batch import image_path, reference_path, fasta_res_group
from hailtop.batch import Batch, ResourceFile
from hailtop.batch.job import Job

from cpg_pipes import Path, to_path
from cpg_pipes.hb.command import wrap_command, seds_to_extend_sample_ids
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.jobs.scripts import check_pedigree
from cpg_pipes.slack import slack_message_cmd, slack_env
from cpg_pipes.targets import Dataset, Sample
from cpg_pipes.types import CramPath, GvcfPath, SequencingType
from cpg_pipes.utils import can_reuse

logger = logging.getLogger(__file__)


# We exclude contaminated samples from relatedness checks
MIN_FREEMIX = 0.04


def pedigree(
    b,
    dataset: Dataset,
    expected_ped_path: Path,
    input_path_by_sid: dict[str, Path | str],
    overwrite: bool,
    out_samples_path: Path,
    out_pairs_path: Path,
    out_html_path: Path,
    out_html_url: str,
    out_checks_path: Path | None = None,
    verifybamid_by_sid: dict[str, Path | str] | None = None,
    label: str | None = None,
    ignore_missing: bool = False,
    job_attrs: dict | None = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
    send_to_slack: bool = True,
) -> list[Job]:
    """
    Add somalier and peddy based jobs that infer relatedness and sex, compare that
    to the provided PED file, and attempt to recover it. If unable to recover, cancel
    the further workflow jobs.

    Returns a job, a path to a fixed PED file if able to recover, and a path to a file
    with relatedness information for each sample pair

    input_path_by_sid can have paths to CRAMs, BAMs, GVCFs, or .somalier
    fingerprints, which will be determined based on extention. Unless a .somalier
    print is provided, `somalier extract` will be run to extract one.
    """
    extract_jobs, somalier_file_by_sample = _prep_somalier_files(
        b=b,
        samples=dataset.get_samples(),
        input_path_by_sid=input_path_by_sid,
        overwrite=overwrite,
        label=label,
        ignore_missing=ignore_missing,
        job_attrs=job_attrs,
        sequencing_type=sequencing_type,
    )

    relate_j = _relate(
        b=b,
        somalier_file_by_sid=somalier_file_by_sample,
        verifybamid_by_sid=verifybamid_by_sid,
        sample_ids=dataset.get_sample_ids(),
        rich_id_map=dataset.rich_id_map(),
        expected_ped_path=expected_ped_path,
        label=label,
        extract_jobs=extract_jobs,
        out_samples_path=out_samples_path,
        out_pairs_path=out_pairs_path,
        out_html_path=out_html_path,
        out_html_url=out_html_url,
        job_attrs=job_attrs,
    )

    check_j = check_pedigree_job(
        b=b,
        samples_file=relate_j.output_samples,
        pairs_file=relate_j.output_pairs,
        expected_ped=b.read_input(str(expected_ped_path)),
        somalier_html_url=out_html_url,
        rich_id_map=dataset.rich_id_map(),
        dataset_name=dataset.name,
        label=label,
        out_checks_path=out_checks_path,
        job_attrs=job_attrs,
        send_to_slack=send_to_slack,
    )
    check_j.depends_on(relate_j)

    return extract_jobs + [relate_j, check_j]


def _make_sample_map(dataset: Dataset):
    """
    Creating sample map to remap internal IDs to participant IDs
    """
    sample_map_fpath = dataset.tmp_prefix() / 'pedigree' / 'sample_map.tsv'
    df = pd.DataFrame(
        [{'id': s.id, 'pid': s.participant_id} for s in dataset.get_samples()]
    )
    with sample_map_fpath.open('w') as fp:
        df.to_csv(fp, sep='\t', index=False, header=False)
    return sample_map_fpath


def ancestry(
    b,
    dataset: Dataset,
    input_path_by_sid: dict[str, Path | str],
    out_tsv_path: Path,
    out_html_path: Path,
    overwrite: bool = True,
    out_html_url: str | None = None,
    label: str | None = None,
    ignore_missing: bool = False,
    job_attrs: dict | None = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
) -> Job:
    """
    Run somalier ancestry https://github.com/brentp/somalier/wiki/ancestry
    """
    extract_jobs, somalier_file_by_sample = _prep_somalier_files(
        b=b,
        samples=dataset.get_samples(),
        input_path_by_sid=input_path_by_sid,
        overwrite=overwrite,
        label=label,
        ignore_missing=ignore_missing,
        job_attrs=job_attrs,
        sequencing_type=sequencing_type,
    )
    j = _ancestry(
        b=b,
        somalier_file_by_sample=somalier_file_by_sample,
        sample_ids=dataset.get_sample_ids(),
        external_id_map=dataset.rich_id_map(),
        label=label,
        extract_jobs=extract_jobs,
        out_tsv_path=out_tsv_path,
        out_html_path=out_html_path,
        out_html_url=out_html_url,
        job_attrs=job_attrs,
    )
    if out_html_url:
        slack_message_cmd(
            j,
            text=(
                f'*[{dataset.name}]* ancestry report: '
                f'<{out_html_url}|{basename(out_html_url)}>'
            ),
        )
    return j


def _prep_somalier_files(
    b,
    samples: list[Sample],
    input_path_by_sid: dict[str, Path | str],
    overwrite: bool,
    label: str | None = None,
    ignore_missing: bool = False,
    job_attrs: dict | None = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
) -> tuple[list[Job], dict[str, Path]]:
    """
    Generate .somalier file for each input.
    """
    extract_jobs = []
    missing_input = []
    somalier_file_by_sample = dict()
    for sample in samples:
        input_path = input_path_by_sid.get(sample.id)
        if input_path is None:
            missing_input.append(sample)
            logger.error(f'Not found somalier input for {sample.id}')
            continue

        input_path = to_path(input_path)

        if input_path.name.endswith('.somalier'):
            somalier_file_by_sample[sample.id] = input_path
            continue

        gvcf_or_cram_or_bam_path: CramPath | GvcfPath
        if input_path.name.endswith('.cram') or input_path.name.endswith('.bam'):
            gvcf_or_cram_or_bam_path = CramPath(input_path)
        else:
            gvcf_or_cram_or_bam_path = GvcfPath(input_path)
        j = extact_job(
            b=b,
            gvcf_or_cram_or_bam_path=gvcf_or_cram_or_bam_path,
            overwrite=overwrite,
            label=label,
            output_path=gvcf_or_cram_or_bam_path.somalier_path,
            job_attrs=(job_attrs or {}) | sample.get_job_attrs(),
            sequencing_type=sequencing_type,
        )
        somalier_file_by_sample[sample.id] = gvcf_or_cram_or_bam_path.somalier_path
        extract_jobs.append(j)

    if len(missing_input) > 0:
        msg = (
            f'Could not find input for ' f'{len(missing_input)}/{len(samples)} samples'
        )
        if ignore_missing:
            logger.warning(msg)
        else:
            raise ValueError(msg)
    return extract_jobs, somalier_file_by_sample


def check_pedigree_job(
    b: Batch,
    samples_file: ResourceFile,
    pairs_file: ResourceFile,
    expected_ped: ResourceFile,
    somalier_html_url: str,
    dataset_name: str,
    rich_id_map: dict[str, str] | None = None,
    label: str | None = None,
    out_checks_path: Path | None = None,
    job_attrs: dict | None = None,
    send_to_slack: bool = True,
) -> Job:
    """
    Run job that checks pedigree and batch correctness. The job will send a Slack
    message about any mismatches.
    """
    title = 'Pedigree check'
    if label:
        title += f' [{label}]'

    check_j = b.new_job(title, job_attrs)
    check_j.attributes['tool'] = 'python'
    STANDARD.set_resources(check_j, ncpu=2)
    check_j.image(image_path('peddy'))

    script_path = to_path(check_pedigree.__file__)
    script_name = script_path.name
    cmd = f"""\
    {seds_to_extend_sample_ids(rich_id_map, [samples_file, pairs_file, expected_ped])
    if rich_id_map else ''}
    python3 {script_name} \\
    --somalier-samples {samples_file} \\
    --somalier-pairs {pairs_file} \\
    --expected-ped {expected_ped} \\
    --html-url {somalier_html_url} \\
    --dataset {dataset_name} \\
    --title "{title}" \\
    --{"no-" if not send_to_slack else ""}send-to-slack

    touch {check_j.output}
    echo "HTML URL: {somalier_html_url}"
    """

    slack_env(check_j)
    check_j.command(
        wrap_command(
            cmd,
            python_script_path=script_path,
            setup_gcp=True,
        )
    )
    if out_checks_path:
        b.write_output(check_j.output, str(out_checks_path))
    return check_j


def _ancestry(
    b: Batch,
    somalier_file_by_sample: dict[str, Path],
    sample_ids: list[str],
    external_id_map: dict[str, str],
    label: str | None,
    extract_jobs: list[Job],
    out_tsv_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    job_attrs: dict | None = None,
) -> Job:
    j = b.new_job('Somalier ancestry' + (f' {label}' if label else ''), job_attrs)
    j.image(image_path('somalier'))
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    STANDARD.set_resources(j, storage_gb=1 + len(sample_ids) // 4000 * 1)
    j.depends_on(*extract_jobs)

    cmd = f"""\
    mkdir /io/batch/1kg
    mv {b.read_input(str(reference_path('somalier_1kg_targz')))} /io/batch/1kg
    (cd /io/batch/1kg && tar -xzf *.tar.gz)

    mkdir /io/batch/somaliers
    """

    for sample_id in sample_ids:
        somalier_file = b.read_input(str(somalier_file_by_sample[sample_id]))
        cmd += f'    cp {somalier_file} /io/batch/somaliers/\n'

    cmd += f"""\
    somalier ancestry \\
    --labels {b.read_input(str(reference_path('somalier_1kg_labels')))} \\
    /io/batch/1kg/1kg-somalier/*.somalier ++ \\
    /io/batch/somaliers/*.somalier \\
    -o ancestry
    ls
    mv ancestry.somalier-ancestry.tsv {j.output_tsv}
    {seds_to_extend_sample_ids(external_id_map, ['ancestry.somalier-ancestry.html'])}
    mv ancestry.somalier-ancestry.html {j.output_html}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'
    j.command(wrap_command(cmd))
    b.write_output(j.output_tsv, str(out_tsv_path))
    b.write_output(j.output_html, str(out_html_path))
    return j


def _relate(
    b: Batch,
    somalier_file_by_sid: dict[str, Path],
    sample_ids: list[str],
    rich_id_map: dict[str, str],
    expected_ped_path: Path,
    label: str | None,
    extract_jobs: list[Job],
    out_samples_path: Path,
    out_pairs_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    verifybamid_by_sid: dict[str, Path] | None = None,
    job_attrs: dict | None = None,
) -> Job:
    title = 'Somalier relate'
    if label:
        title += f' [{label}]'

    j = b.new_job(title, job_attrs)
    j.attributes['tool'] = 'somalier'
    j.image(image_path('somalier'))
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    STANDARD.set_resources(j, storage_gb=1 + len(sample_ids) // 4000 * 1)
    j.depends_on(*extract_jobs)

    cmd = ''
    input_files_file = '/io/input_files.list'
    samples_ids_file = '/io/sample_ids.list'
    cmd += f'touch {input_files_file}'
    cmd += f'touch {samples_ids_file}'
    for sample_id in sample_ids:
        if verifybamid_by_sid:
            if sample_id not in verifybamid_by_sid:
                continue
            somalier_file = b.read_input(str(somalier_file_by_sid[sample_id]))
            cmd += f"""
            FREEMIX=$(cat {b.read_input(str(verifybamid_by_sid[sample_id]))} | tail -n1 | cut -f7)
            if [[ $(echo "$FREEMIX > {MIN_FREEMIX}" | bc) -eq 0 ]]; \
            then echo "{somalier_file}" >> {input_files_file}; \
            echo "{sample_id}" >> {samples_ids_file}; fi
            """
        else:
            somalier_file = b.read_input(str(somalier_file_by_sid[sample_id]))
            cmd += f"""
            echo "{somalier_file}" >> {input_files_file}
            echo "{sample_id}" >> {samples_ids_file}
            """

    cmd += f"""
    cat {b.read_input(str(expected_ped_path))} | \
    grep -v Family.ID | \
    grep -f {samples_ids_file} > expected.ped 
    """

    cmd += f"""
    somalier relate \\
    $(cat {input_files_file}) \\
    --ped expected.ped \\
    -o related \\
    --infer
    ls
    mv related.pairs.tsv {j.output_pairs}
    mv related.samples.tsv {j.output_samples}
    {seds_to_extend_sample_ids(rich_id_map, ['related.html'])}
    mv related.html {j.output_html}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'

    j.command(wrap_command(cmd))
    # Copy somalier outputs to final destination.
    b.write_output(j.output_samples, str(out_samples_path))
    b.write_output(j.output_pairs, str(out_pairs_path))
    b.write_output(j.output_html, str(out_html_path))
    return j


def extact_job(
    b,
    gvcf_or_cram_or_bam_path: CramPath | GvcfPath,
    overwrite: bool = True,
    label: str | None = None,
    output_path: Path | None = None,
    job_attrs: dict | None = None,
    sequencing_type: SequencingType = SequencingType.GENOME,
) -> Job | None:
    """
    Run "somalier extract" to generate a fingerprint for a `sample`
    from `fpath` (which can be a GVCF, a CRAM or a BAM)
    """
    if can_reuse(output_path, overwrite):
        return None

    job_attrs = (job_attrs or {}) | {'tool': 'somalier'}
    j = b.new_job('Somalier extract' + (f' {label}' if label else ''), job_attrs)

    if not output_path:
        output_path = gvcf_or_cram_or_bam_path.somalier_path

    j.image(image_path('somalier'))
    if isinstance(gvcf_or_cram_or_bam_path, CramPath):
        storage_gb = None  # avoid extra disk by default
        if sequencing_type == SequencingType.GENOME:
            storage_gb = 100
            if gvcf_or_cram_or_bam_path.is_bam:
                storage_gb = 200
        STANDARD.set_resources(j, ncpu=4, storage_gb=storage_gb)
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.index_path),
        )
    else:
        STANDARD.set_resources(j, ncpu=2, storage_gb=10)
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.tbi_path),
        )

    ref = fasta_res_group(b)
    sites = b.read_input(str(reference_path('somalier_sites')))

    cmd = f"""\
    SITES=/io/batch/sites/{reference_path('somalier_sites').name}
    retry gsutil cp {reference_path('somalier_sites')} $SITES

    somalier extract -d extracted/ --sites {sites} -f {ref.base} \\
    {input_file['base']}

    mv extracted/*.somalier {j.output_file}
    """
    j.command(wrap_command(cmd, setup_gcp=True, define_retry_function=True))
    b.write_output(j.output_file, str(output_path))
    return j
