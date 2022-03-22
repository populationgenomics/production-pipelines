"""
Adding jobs for fingerprinting and pedigree checks. Mostly using Somalier.
"""
import logging
import pandas as pd
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import Path, to_path
from cpg_pipes import images, utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.types import CramPath, GvcfPath
from cpg_pipes.pipeline.targets import Dataset, Sample
from cpg_pipes.refdata import RefData

logger = logging.getLogger(__file__)


def pedigree(
    b,
    dataset: Dataset,
    input_path_by_sid: dict[str, Path | str],
    refs: RefData,
    overwrite: bool,
    out_samples_path: Path | None = None,
    out_pairs_path: Path | None = None,
    out_html_path: Path | None = None,
    out_html_url: str | None = None,
    label: str | None = None,
    ignore_missing: bool = False,
    dry_run: bool = False,
) -> Job:
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
        dataset=dataset,
        refs=refs,
        input_path_by_sid=input_path_by_sid,
        overwrite=overwrite,
        label=label,
        ignore_missing=ignore_missing,
    )
    
    relate_j = _relate(
        b=b,
        somalier_file_by_sample=somalier_file_by_sample,
        dataset=dataset,
        label=label,
        extract_jobs=extract_jobs,
        out_samples_path=out_samples_path,
        out_pairs_path=out_pairs_path,
        out_html_path=out_html_path,
        out_html_url=out_html_url,
        dry_run=dry_run,
    )
    
    _check_pedigree(
        b=b,
        dataset=dataset,
        relate_j=relate_j,
        label=label,
        somalier_html_url=out_html_url,
        dry_run=dry_run,
    )

    return relate_j


def ancestry(
    b,
    dataset: Dataset,
    refs: RefData,
    input_path_by_sid: dict[str, Path | str],
    overwrite: bool,
    out_tsv_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
    label: str | None = None,
    ignore_missing: bool = False,
) -> Job:
    """
    Run somalier ancestry https://github.com/brentp/somalier/wiki/ancestry
    """
    extract_jobs, somalier_file_by_sample = _prep_somalier_files(
        b=b,
        dataset=dataset,
        refs=refs,
        input_path_by_sid=input_path_by_sid,
        overwrite=overwrite,
        label=label,
        ignore_missing=ignore_missing,
    )

    j = _ancestry(
        b=b,
        somalier_file_by_sample=somalier_file_by_sample,
        dataset=dataset,
        refs=refs,
        label=label,
        extract_jobs=extract_jobs,
        out_tsv_path=out_tsv_path,
        out_html_path=out_html_path,
        out_html_url=out_html_url,
    )
    return j


def _prep_somalier_files(
    b,
    dataset: Dataset,
    input_path_by_sid: dict[str, Path | str],
    refs: RefData,
    overwrite: bool,
    label: str | None = None,
    ignore_missing: bool = False,
) -> tuple[list[Job], dict[str, Path]]:
    """
    Generate .somalier file for each input
    """
    extract_jobs = []
    missing_input = []
    somalier_file_by_sample = dict()
    for sample in dataset.get_samples():
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
            sample=sample,
            refs=refs,
            gvcf_or_cram_or_bam_path=gvcf_or_cram_or_bam_path,
            overwrite=overwrite,
            label=label,
            out_fpath=gvcf_or_cram_or_bam_path.somalier_path,
        )
        somalier_file_by_sample[sample.id] = gvcf_or_cram_or_bam_path.somalier_path
        extract_jobs.append(j)

    if len(missing_input) > 0:
        msg = (
            f'Could not find input for '
            f'{len(missing_input)}/{len(dataset.get_samples())} samples'
        )
        if ignore_missing:
            logger.warning(msg)
        else:
            raise ValueError(msg)
    return extract_jobs, somalier_file_by_sample


def _check_pedigree(
    b: Batch, 
    dataset: Dataset, 
    relate_j: Job, 
    label: str | None,
    somalier_html_url: str | None = None,
    dry_run: bool = False,
):
    check_j = b.new_job(
        'Check relatedness and sex' + (f' {label}' if label else ''),
        dict(dataset=dataset.name),
    )
    STANDARD.set_resources(check_j, ncpu=2)
    check_j.image(images.PEDDY_IMAGE)
    # Creating sample map to remap internal IDs to participant IDs
    sample_map_fpath = dataset.get_tmp_bucket() / 'pedigree' / 'sample_maps' / f'{dataset.name}.tsv'
    if not dry_run:
        df = pd.DataFrame([
            {'id': s.id, 'pid': s.participant_id} for s in dataset.get_samples()
        ])
        df.to_csv(str(sample_map_fpath), sep='\t', index=False, header=False)
        
    script_name = 'check_pedigree.py'
    script_path = to_path(__file__).parent.parent.parent / utils.SCRIPTS_DIR / script_name
    with open(script_path) as f:
        script = f.read()
    # We do not wrap the command nicely to avoid breaking python indents of {script}
    cmd = f"""\
    cat <<EOT >> {script_name}
    {script}
    EOT
    python {script_name} \
    --somalier-samples {relate_j.output_samples} \
    --somalier-pairs {relate_j.output_pairs} \
    --sample-map {b.read_input(str(sample_map_fpath))} \
    {('--somalier-html ' + somalier_html_url) if somalier_html_url else ''}
    """
    check_j.command(cmd)
    check_j.depends_on(relate_j)


def _ancestry(
    b: Batch, 
    somalier_file_by_sample: dict[str, Path],
    dataset: Dataset,
    refs: RefData,
    label: str | None,
    extract_jobs: list[Job],
    out_tsv_path: Path,
    out_html_path: Path,
    out_html_url: str | None = None,
) -> Job:
    j = b.new_job(
        'Somalier ancestry' + (f' {label}' if label else ''),
        dict(dataset=dataset.name),
    )
    j.image(images.BIOINFO_IMAGE)
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    STANDARD.set_resources(
        j, 
        storage_gb=1 + len(dataset.get_samples()) // 4000 * 1,
    )
    j.depends_on(*extract_jobs)

    cmd = f"""\
    mkdir /io/batch/1kg
    mv {b.read_input(refs.somalier_1kg_targz)} /io/batch/1kg
    (cd /io/batch/1kg && tar -xzf *.tar.gz)

    mkdir /io/batch/somaliers
    """

    for sample in dataset.get_samples():
        if sample.id:
            somalier_file = b.read_input(str(somalier_file_by_sample[sample.id]))
            cmd += f'    cp {somalier_file} /io/batch/somaliers/\n'

    cmd += f"""\
    somalier ancestry \\
    --labels {b.read_input(refs.somalier_1kg_labels_tsv)} \\
    /io/batch/1kg/1kg-somalier/*.somalier ++ \\
    /io/batch/somaliers/*.somalier \\
    -o ancestry

    mv ancestry.somalier-ancestry.tsv {j.output_tsv}
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
    somalier_file_by_sample: dict[str, Path],
    dataset: Dataset,
    label: str | None,
    extract_jobs: list[Job],
    out_samples_path: Path | None = None,
    out_pairs_path: Path | None = None,
    out_html_path: Path | None = None,
    out_html_url: str | None = None,
    dry_run: bool = False,
) -> Job:
    j = b.new_job(
        'Somalier relate' + (f' {label}' if label else ''),
        dict(dataset=dataset.name),
    )
    j.image(images.BIOINFO_IMAGE)
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    STANDARD.set_resources(
        j, 
        storage_gb=1 + len(dataset.get_samples()) // 4000 * 1,
    )
    j.depends_on(*extract_jobs)
    ped_fpath = dataset.get_tmp_bucket() / f'{dataset.name}.ped'
    datas = []
    for sample in dataset.get_samples():
        if sample.pedigree:
            datas.append(sample.pedigree.get_ped_dict())
    df = pd.DataFrame(datas)
    if not dry_run:
        df.to_csv(str(ped_fpath), sep='\t', index=False)
    ped_file = b.read_input(str(ped_fpath))
    input_files_lines = ''
    for sample in dataset.get_samples():
        if sample.id:
            somalier_file = b.read_input(str(somalier_file_by_sample[sample.id]))
            input_files_lines += f'{somalier_file} \\\n'
    cmd = f"""\
    cat {ped_file} | grep -v Family.ID > /io/samples.ped 
    
    somalier relate \\
    {input_files_lines} \\
    --ped /io/samples.ped \\
    -o related \\
    --infer
    
    ls
    mv related.html {j.output_html}
    mv related.pairs.tsv {j.output_pairs}
    mv related.samples.tsv {j.output_samples}
    """
    if out_html_url:
        cmd += '\n' + f'echo "HTML URL: {out_html_url}"'
    j.command(wrap_command(cmd))
    # Copy somalier outputs to buckets
    b.write_output(j.output_samples, str(out_samples_path))
    b.write_output(j.output_pairs, str(out_pairs_path))
    if out_html_path:
        b.write_output(j.output_html, str(out_html_path))
    return j


def extact_job(
    b,
    sample: Sample,
    gvcf_or_cram_or_bam_path: CramPath | GvcfPath,
    refs: RefData,
    overwrite: bool,
    label: str | None = None,
    out_fpath: Path | None = None,
) -> Job:
    """
    Run "somalier extract" to generate a fingerprint for a `sample`
    from `fpath` (which can be a GVCF, a CRAM or a BAM)
    """
    j = b.new_job(
        'Somalier extract' + (f' {label}' if label else ''),
        sample.get_job_attrs(),
    )

    if not out_fpath:
        out_fpath = gvcf_or_cram_or_bam_path.somalier_path

    if utils.can_reuse(out_fpath, overwrite):
        j.name += ' [reuse]'
        return j

    j.image(images.BIOINFO_IMAGE)
    j.memory('standard')
    if isinstance(gvcf_or_cram_or_bam_path, CramPath):
        STANDARD.request_resources(
            ncpu=4, 
            storage_gb=200 if gvcf_or_cram_or_bam_path.is_bam else 50
        )
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.index_path)
        )
    else:
        STANDARD.request_resources(ncpu=2, storage_gb=10)
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.tbi_path)
        )
    
    ref = refs.fasta_res_group(b)
    sites = b.read_input(refs.somalier_sites)

    cmd = f"""\
    SITES=/io/batch/sites/{refs.somalier_sites.name}
    retry gsutil cp {refs.somalier_sites} $SITES

    somalier extract -d extracted/ --sites {sites} -f {ref.base} \\
    {input_file['base']}

    mv extracted/*.somalier {j.output_file}
    """
    j.command(wrap_command(cmd, setup_gcp=True, define_retry_function=True))
    b.write_output(j.output_file, str(out_fpath))
    return j
