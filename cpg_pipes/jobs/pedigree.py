import sys
from os.path import join, dirname, pardir, basename
from typing import Optional, List, Tuple, Dict
import logging

from hailtop.batch.job import Job
import pandas as pd

from cpg_pipes import utils, resources
from cpg_pipes.pipeline import Project, Sample
from cpg_pipes.hailbatch import wrap_command

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


def add_pedigree_jobs(
    b,
    project: Project,
    input_path_by_sid: Dict[str, str],
    overwrite: bool,
    fingerprints_bucket: str,
    tmp_bucket: str,
    local_tmp_dir: str,
    web_bucket: Optional[str] = None,
    web_url: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    label: str = None,
    ignore_missing: bool = False,
) -> Tuple[Job, str, str]:
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
    for sample in project.get_samples():
        input_path = input_path_by_sid.get(sample.id)
        if input_path is None:
            missing_input.append(sample)
            logger.error(f'Not found somalier input for {sample.id}')
            continue
        
        if input_path.endswith('.somalier'):
            somalier_file_by_sample[sample.id] = input_path
            continue

        j, out_fpath = somalier_extact_job(
            b=b,
            sample=sample,
            gvcf_or_cram_or_bam_path=input_path,
            overwrite=overwrite,
            label=label,
            depends_on=depends_on,
        )
        somalier_file_by_sample[sample.id] = out_fpath
        extract_jobs.append(j)

    if len(missing_input) > 0:
        (logger.critical if not ignore_missing else logger.warning)(
            f'Could not find input for '
            f'{len(missing_input)}/{len(project.get_samples())} samples'
        )
        if not ignore_missing:
            sys.exit(1)

    relate_j = b.new_job(
        'Somalier relate' + (f' {label}' if label else ''), 
        dict(project=project.name),
    )
    relate_j.image(resources.BIOINFO_IMAGE)
    relate_j.cpu(1)
    relate_j.memory('standard')  # ~ 4G/core ~ 4G
    # Size of one somalier file is 212K, so we add another G only if the number of
    # samples is >4k
    relate_j.storage(f'{1 + len(project.get_samples()) // 4000 * 1}G')
    if depends_on:
        extract_jobs.extend(depends_on)
    relate_j.depends_on(*extract_jobs)

    ped_fpath = join(tmp_bucket, f'{project.name}.ped')
    datas = []
    for sample in project.get_samples():
        if sample.pedigree:
            datas.append(sample.pedigree.get_ped_dict())
    df = pd.DataFrame(datas)
    df.to_csv(ped_fpath, sep='\t', index=False)
    ped_file = b.read_input(ped_fpath)

    input_files_lines = ''
    for sample in project.get_samples():
        if sample.id:
            somalier_file = b.read_input(somalier_file_by_sample[sample.id])
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

    # Copy somalier outputs to buckets
    prefix = join(fingerprints_bucket, project.name)
    somalier_samples_path = f'{prefix}.samples.tsv'
    somalier_pairs_path = f'{prefix}.pairs.tsv'
    b.write_output(relate_j.output_samples, somalier_samples_path)
    b.write_output(relate_j.output_pairs, somalier_pairs_path)
    # Copy somalier HTML to the web bucket
    somalier_html_url = None
    if web_bucket and web_url:
        rel_path = f'pedigree/{project.name}.html'
        somalier_html_path = f'{web_bucket}/{rel_path}'
        somalier_html_url = f'{web_url}/{rel_path}'
        b.write_output(relate_j.output_html, somalier_html_path)

    check_j = b.new_job(
        'Check relatedness and sex' + (f' {label}' if label else ''), 
        dict(project=project.name),
    )
    check_j.image(resources.PEDDY_IMAGE)
    check_j.cpu(1)
    check_j.memory('standard')  # ~ 3.75G/core ~ 3.75G

    local_sample_map_fpath = join(local_tmp_dir, f'{project.name}-sample-map.tsv')
    remote_sample_map_fpath = join(tmp_bucket, f'{project.name}-sample-map.tsv')
    with open(local_sample_map_fpath, 'w') as f:
        for sample in project.get_samples():
            f.write(f'{sample.id}\t{sample.participant_id}\n')
    utils.gsutil_cp(local_sample_map_fpath, remote_sample_map_fpath)

    script_name = 'check_pedigree.py'
    script_path = join(dirname(__file__), pardir, pardir, utils.SCRIPTS_DIR, script_name)
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
--sample-map {b.read_input(remote_sample_map_fpath)} \
{('--somalier-html ' + somalier_html_url) if somalier_html_url else ''}
"""
    check_j.command(cmd)

    check_j.depends_on(relate_j)
    return relate_j, somalier_samples_path, somalier_pairs_path


def somalier_extact_job(
    b,
    sample: Sample,
    gvcf_or_cram_or_bam_path: str,
    overwrite: bool,
    label: Optional[str] = None,
    depends_on: Optional[List[Job]] = None,
    out_fpath: Optional[str] = None,
) -> Tuple[Job, str]:
    """
    Run "somalier extract" to generate a fingerprint for a `sample`
    from `fpath` (which can be a gvcf, a cram or a bam)
    """
    j = b.new_job(
        'Somalier extract' + (f' {label}' if label else ''),
        dict(sample=sample.id, project=sample.project.name),
    )

    if not out_fpath:
        out_fpath = gvcf_or_cram_or_bam_path\
            .replace('.cram', '.somalier')\
            .replace('.bam', '.somalier')\
            .replace('.g.vcf.gz', '.somalier')

    if utils.can_reuse(out_fpath, overwrite):
        j.name += ' [reuse]'
        return j, out_fpath

    j.image(resources.BIOINFO_IMAGE)
    j.memory('standard')
    if gvcf_or_cram_or_bam_path.endswith('.bam'):
        j.cpu(4)
        j.storage(f'200G')
        input_file = b.read_input_group(
            base=gvcf_or_cram_or_bam_path,
            index=gvcf_or_cram_or_bam_path + '.bai',
        )
    elif gvcf_or_cram_or_bam_path.endswith('.cram'):
        j.cpu(4)
        j.storage(f'50G')
        input_file = b.read_input_group(
            base=gvcf_or_cram_or_bam_path,
            index=gvcf_or_cram_or_bam_path + '.crai',
        )
    else:
        j.cpu(2)
        j.storage(f'10G')
        input_file = b.read_input_group(
            base=gvcf_or_cram_or_bam_path,
            index=gvcf_or_cram_or_bam_path + '.tbi',
        )
    
    if depends_on:
        j.depends_on(*depends_on)

    ref_fasta = resources.REF_FASTA
    ref_fai = resources.REF_FASTA + '.fai'
    ref_dict = (
        ref_fasta.replace('.fasta', '').replace('.fna', '').replace('.fa', '') + '.dict'
    )

    cmd = f"""\
    # Copying reference data to avoid GCP bandwidth limits
    retry gsutil cp {ref_fasta} /io/batch/{basename(ref_fasta)}
    retry gsutil cp {ref_fai}   /io/batch/{basename(ref_fai)}
    retry gsutil cp {ref_dict}  /io/batch/{basename(ref_dict)}
    SITES=/io/batch/sites/{basename(resources.SOMALIER_SITES)}
    retry gsutil cp {resources.SOMALIER_SITES} $SITES

    somalier extract -d extracted/ --sites $SITES \\
    -f /io/batch/{basename(ref_fasta)} \\
    {input_file['base']}

    mv extracted/*.somalier {j.output_file}
    """
    j.command(wrap_command(cmd, setup_gcp=True, define_retry_function=True))
    b.write_output(j.output_file, out_fpath)
    return j, out_fpath
