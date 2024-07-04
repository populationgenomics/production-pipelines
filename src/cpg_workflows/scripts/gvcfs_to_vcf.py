#! /usr/bin/env python3

"""
all the steps from gvcfs to vcf
required gcloud authentication
"""


from argparse import ArgumentParser
from os.path import join
from subprocess import run

import hail as hl

from cpg_utils.config import get_config
from cpg_utils.hail_batch import get_batch, image_path, init_batch, output_path
from cpg_workflows.utils import chunks, get_logger

MAX_COMPOSE: int = 30


def gvcfs_to_vds(gvcfs: list[str], vds_out: str, sgids: list[str], family: str):
    """
    start up a combiner instance, and combine the gVCFs

    Args:
        gvcfs (list[str]):
        vds_out (str):
        sgids (list[str]):
        family (str):
    """

    tmp_prefix = output_path(f'{family}.vds', category='tmp')

    if get_config()['workflow']['sequencing_type'] == 'exome':
        params = {'use_exome_default_intervals': True}
    else:
        params = {'use_genome_default_intervals': True}

    combiner = hl.vds.new_combiner(
        gvcf_paths=gvcfs,
        gvcf_sample_names=sgids,
        gvcf_external_header=gvcfs[0],
        output_path=vds_out,
        reference_genome='GRCh38',
        temp_path=tmp_prefix,
        **params,
    )
    combiner.run()


def vds_to_vcf(vds_path: str, output: str):
    vds = hl.vds.read_vds(vds_path)
    vds = hl.vds.split_multi(vds)
    mt = hl.vds.to_dense_mt(vds)

    # gotta drop this (it's a dict)
    if 'gvcf_info' in mt.entry:
        mt = mt.drop('gvcf_info')

    hl.export_vcf(mt, output, parallel='separate_header')
    get_logger().info(f'Exported VCF fragments to {output}')


def squash_fragments_to_vcf(vcf_fragment_dir: str, vcf_out: str):
    """
    Squash VCF fragments into a single VCF file

    Args:
        vcf_fragment_dir ():
        vcf_out ():

    Returns:

    """

    get_logger().info('Parsing the manifest file and generating a bash script for concatenating the VCF files')
    manifest = hl.hadoop_open(join(vcf_fragment_dir, 'shard-manifest.txt')).readlines()

    # prefix these shard names to get full GCP path for each
    manifest_paths = [join(vcf_fragment_dir, str(fragment).strip()) for fragment in manifest]

    # maybe don't do this
    if len(manifest) <= MAX_COMPOSE:
        run(['gcloud', 'storage', 'objects', 'compose', *manifest_paths, vcf_out])
        get_logger().info(f'Composed {len(manifest)} VCF fragments into {vcf_out}')
        return

    # rolling squash of the chunks, should enable infinite-ish scaling
    temp_chunk_prefix_num = 1
    while len(manifest_paths) > 1:

        # ordering is super important
        chunks_this_round = []
        condense_temp = join(vcf_fragment_dir, f'temp_chunk_{temp_chunk_prefix_num}')
        temp_chunk_prefix_num += 1

        for i, sub_chunk in enumerate(chunks(manifest_paths, MAX_COMPOSE)):
            sub_chunk_out = join(condense_temp, f'chunk_{i}.vcf.bgz')
            run(['gcloud', 'storage', 'objects', 'compose', *sub_chunk, sub_chunk_out])
            chunks_this_round.append(sub_chunk_out)

        manifest_paths = chunks_this_round

    # tabix the results - read the final vcf in, index it, move the index, and write it out
    input_vcf = get_batch().read_input(manifest_paths[0])
    final_job = get_batch().new_bash_job(name='index_final_vcf')
    final_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
    final_job.image(image_path('bcftools'))
    final_job.storage('10Gi')  # make this configurable
    final_job.command(f'mv {input_vcf} {final_job.output["vcf.bgz"]} && tabix {final_job.output["vcf.bgz"]}')
    get_batch().write_output(final_job.output, vcf_out.removesuffix('.vcf.bgz'))

    # this will be a batch-in-batch, and should execute quickly
    get_batch().run(wait=True)


if __name__ == '__main__':

    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--gvcfs', help='Input gVCFs', nargs='+')
    parser.add_argument('--sgids', help='Family Member IDs', nargs='+')
    parser.add_argument('--out', help='Output VCF path')
    parser.add_argument('--family', help='Family Name/Number')
    args = parser.parse_args()

    assert len(args.gvcfs) == len(args.sgids)

    init_batch()

    vds_path = output_path(f'{args.family}.vds', category='tmp')

    get_logger(__file__).info(f'Creating VCF {args.out} from {len(args.gvcfs)} gVCFs')

    gvcfs_to_vds(gvcfs=args.gvcfs, sgids=args.sgids, vds_out=vds_path, family=args.family)

    get_logger(__file__).info('Creating VCF fragments from VDS...')

    vcf_fragments_tmp = output_path(f'fragments_{args.family}.vcf.bgz', category='tmp')

    vds_to_vcf(vds_path=vds_path, output=vcf_fragments_tmp)

    get_logger(__file__).info('Creating single VCF from fragments...')

    squash_fragments_to_vcf(vcf_fragment_dir=vcf_fragments_tmp, vcf_out=args.out)

    get_logger(__file__).info(f'All done writing {args.out}...')
