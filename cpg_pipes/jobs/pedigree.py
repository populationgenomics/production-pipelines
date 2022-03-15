"""
Adding jobs for fingerprinting and pedigree checks. Mostly using Somalier.
"""
import os
import sys
import logging
from pathlib import Path
from cloudpathlib import CloudPath
import pandas as pd
from hailtop.batch.job import Job
from hailtop.batch import Batch

from cpg_pipes import buckets, images, ref_data, utils
from cpg_pipes.hb.command import wrap_command
from cpg_pipes.hb.resources import STANDARD
from cpg_pipes.pipeline.analysis import CramPath, GvcfPath
from cpg_pipes.pipeline.dataset import Dataset
from cpg_pipes.pipeline.sample import Sample

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def add_pedigree_jobs(
    b,
    dataset: Dataset,
    input_path_by_sid: dict[str, CloudPath | str],
    overwrite: bool,
    tmp_bucket: CloudPath,
    out_samples_path: CloudPath | None = None,
    out_pairs_path: CloudPath | None = None,
    out_ancestry_path: CloudPath | None = None,
    out_html_path: CloudPath | None = None,
    fingerprints_bucket: CloudPath | None = None,
    web_bucket: CloudPath | None = None,
    web_url: str | None = None,
    depends_on: list[Job] | None = None,
    label: str | None = None,
    ignore_missing: bool = False,
    dry_run: bool = False,
) -> tuple[Job, CloudPath, CloudPath]:
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
    extract_jobs = []
    missing_input = []
    somalier_file_by_sample = dict()
    for sample in dataset.get_samples():
        input_path = input_path_by_sid.get(sample.id)
        if input_path is None:
            missing_input.append(sample)
            logger.error(f'Not found somalier input for {sample.id}')
            continue
            
        input_path = CloudPath(input_path)

        if input_path.name.endswith('.somalier'):
            somalier_file_by_sample[sample.id] = input_path
            continue

        gvcf_or_cram_or_bam_path: CramPath | GvcfPath
        if input_path.name.endswith('.cram') or input_path.name.endswith('.bam'):
            gvcf_or_cram_or_bam_path = CramPath(input_path)
        else:
            gvcf_or_cram_or_bam_path = GvcfPath(input_path)
        j, out_fpath = somalier_extact_job(
            b=b,
            sample=sample,
            gvcf_or_cram_or_bam_path=gvcf_or_cram_or_bam_path,
            overwrite=overwrite,
            label=label,
            depends_on=depends_on,
        )
        somalier_file_by_sample[sample.id] = out_fpath
        extract_jobs.append(j)

    if len(missing_input) > 0:
        (logger.critical if not ignore_missing else logger.warning)(
            f'Could not find input for '
            f'{len(missing_input)}/{len(dataset.get_samples())} samples'
        )
        if not ignore_missing:
            sys.exit(1)

    relate_j = _relate_job(
        b=b,
        somalier_file_by_sample=somalier_file_by_sample,
        dataset=dataset,
        tmp_bucket=tmp_bucket,
        label=label,
        extract_jobs=extract_jobs,
        depends_on=depends_on,
        dry_run=dry_run,
    )

    pairs_path, samples_path, html_url = _copy_somalier_output(
        b=b, 
        dataset=dataset,
        relate_j=relate_j, 
        out_samples_path=out_samples_path,
        out_pairs_path=out_pairs_path,
        out_html_path=out_html_path,
        out_ancestry_path=out_ancestry_path,
        fingerprints_bucket=fingerprints_bucket, 
        web_bucket=web_bucket, 
        web_url=web_url,
    )

    _check_pedigree_job(
        b=b,
        dataset=dataset,
        relate_j=relate_j,
        tmp_bucket=tmp_bucket,
        label=label,
        dry_run=dry_run,
        somalier_html_url=html_url,
    )

    return relate_j, samples_path, pairs_path


def _copy_somalier_output(
    b: Batch,
    dataset: Dataset, 
    relate_j: Job,
    out_samples_path: CloudPath | None = None,
    out_pairs_path: CloudPath | None = None,
    out_html_path: CloudPath | None = None,
    out_ancestry_path: CloudPath | None = None,
    fingerprints_bucket: CloudPath | None = None,
    web_bucket: CloudPath | None = None,
    web_url: str | None = None,
) -> tuple[CloudPath, CloudPath, str | None]:
    """
    Copy somalier outputs to buckets
    """
    if not (out_samples_path and out_pairs_path):
        if not fingerprints_bucket:
            raise ValueError(
                'Somalier Either both out_samples_path/out_pairs_path, or '
                'fingerprints_bucket should be specified'
            )
        prefix = fingerprints_bucket / dataset.name
        out_samples_path = CloudPath(f'{prefix}.samples.tsv')
        out_pairs_path = CloudPath(f'{prefix}.pairs.tsv')
    b.write_output(relate_j.output_samples, str(out_samples_path))
    b.write_output(relate_j.output_pairs, str(out_pairs_path))
    if out_ancestry_path:
        b.write_output(relate_j.output_ancesry, str(out_ancestry_path))

    # Copy somalier HTML to the web bucket
    out_html_url = None
    if not out_html_path:
        if web_bucket:
            out_html_path = web_bucket / 'pedigree' / f'{dataset.name}.html'
        if web_url:
            out_html_url = f'{web_url}/pedigree/{dataset.name}.html'
    else:
        out_html_url = f'{web_url}/{out_html_path.name}'
    
    if out_html_path:
        b.write_output(relate_j.output_html, str(out_html_path))
        
    return out_pairs_path, out_samples_path, out_html_url


def _check_pedigree_job(
    b: Batch, 
    dataset: Dataset, 
    relate_j: Job, 
    tmp_bucket: CloudPath,
    label: str | None,
    dry_run: bool = False,
    somalier_html_url: str | None = None,
):
    check_j = b.new_job(
        'Check relatedness and sex' + (f' {label}' if label else ''),
        dict(dataset=dataset.name),
    )
    STANDARD.set_resources(check_j, ncpu=2)
    check_j.image(images.PEDDY_IMAGE)
    # Creating sample map to remap internal IDs to participant IDs
    sample_map_fpath = tmp_bucket / 'pedigree' / 'sample_maps' / f'{dataset.name}.tsv'
    if not dry_run:
        df = pd.DataFrame([
            {'id': s.id, 'pid': s.participant_id} for s in dataset.get_samples()
        ])
        df.to_csv(str(sample_map_fpath), sep='\t', index=False, header=False)
    script_name = 'check_pedigree.py'
    script_path = Path(__file__).parent.parent.parent / utils.SCRIPTS_DIR / script_name
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


def _relate_job(
    b: Batch, 
    somalier_file_by_sample: dict[str, CloudPath],
    dataset: Dataset,
    tmp_bucket: CloudPath,
    label: str | None,
    extract_jobs: list[Job],
    depends_on: list[Job] | None = None,
    dry_run: bool = False,
) -> Job:
    relate_j = b.new_job(
        'Somalier relate' + (f' {label}' if label else ''),
        dict(dataset=dataset.name),
    )
    relate_j.image(images.BIOINFO_IMAGE)
    relate_j.cpu(1)
    relate_j.memory('standard')  # ~ 4G/core ~ 4G
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    relate_j.storage(f'{1 + len(dataset.get_samples()) // 4000 * 1}G')
    if depends_on:
        extract_jobs.extend(depends_on)
    relate_j.depends_on(*extract_jobs)
    ped_fpath = tmp_bucket / f'{dataset.name}.ped'
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
    mv related.html {relate_j.output_html}
    mv related.pairs.tsv {relate_j.output_pairs}
    mv related.samples.tsv {relate_j.output_samples}
    """
    relate_j.command(wrap_command(cmd))
    return relate_j


def somalier_extact_job(
    b,
    sample: Sample,
    gvcf_or_cram_or_bam_path: CramPath | GvcfPath,
    overwrite: bool,
    label: str | None = None,
    depends_on: list[Job] | None = None,
    out_fpath: CloudPath | None = None,
) -> tuple[Job, CloudPath]:
    """
    Run "somalier extract" to generate a fingerprint for a `sample`
    from `fpath` (which can be a GVCF, a CRAM or a BAM)
    """
    j = b.new_job(
        'Somalier extract' + (f' {label}' if label else ''),
        dict(sample=sample.id, dataset=sample.dataset.name),
    )

    if not out_fpath:
        out_fpath = gvcf_or_cram_or_bam_path.somalier_path

    if buckets.can_reuse(out_fpath, overwrite):
        j.name += ' [reuse]'
        return j, out_fpath

    j.image(images.BIOINFO_IMAGE)
    j.memory('standard')
    if isinstance(gvcf_or_cram_or_bam_path, CramPath):
        j.cpu(4)
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.index_path)
        )
        if gvcf_or_cram_or_bam_path.is_bam:
            j.storage(f'200G')
        else:
            j.storage(f'50G')
    else:
        j.cpu(2)
        j.storage(f'10G')
        input_file = b.read_input_group(
            base=str(gvcf_or_cram_or_bam_path),
            index=str(gvcf_or_cram_or_bam_path.tbi_path)
        )
    
    if depends_on:
        j.depends_on(*depends_on)

    ref_fasta = ref_data.REF_FASTA
    ref_fai = ref_data.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )

    cmd = f"""\
    # Copying reference data to avoid GCP bandwidth limits
    retry gsutil cp {ref_fasta} /io/batch/{os.path.basename(ref_fasta)}
    retry gsutil cp {ref_fai}   /io/batch/{os.path.basename(ref_fai)}
    retry gsutil cp {ref_dict}  /io/batch/{os.path.basename(ref_dict)}
    SITES=/io/batch/sites/{os.path.basename(ref_data.SOMALIER_SITES)}
    retry gsutil cp {ref_data.SOMALIER_SITES} $SITES

    somalier extract -d extracted/ --sites $SITES \\
    -f /io/batch/{os.path.basename(ref_fasta)} \\
    {input_file['base']}

    mv extracted/*.somalier {j.output_file}
    """
    j.command(wrap_command(cmd, setup_gcp=True, define_retry_function=True))
    b.write_output(j.output_file, str(out_fpath))
    return j, out_fpath
