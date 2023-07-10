"""
quick script to find samples we want to analyse

2 arguments:
    - min number of samples to include
    - output path of resulting toml file

e.g.

$ python scripts/gatk_samples_to_run.py 160 samples_here.toml

This can either take a hard-coded list of projects (as below), or query
metamist for 'all projects I have permission to see'.

The script iterates over each project in turn, pulling all valid SG IDs in the
project and grouping by family. It then loops over all the returned data,
selecting a random family and adding all their SGs to a set. Once the set meets
or exceeds the minimum number of samples, the script stops adding new SGs.
If the SGs from this project are exhausted and the minimum requirement has not
been met, it pulls the next project and repeats the process.

When the SGs are being pulled from the metamist results, they are skipped if:
- the SG ID matches a regex that looks like an early/NAGIM sample
- the SG ID is associated with a previously completed GATK-SV analysis
- the SG ID is present in the list of IDs to avoid

Assumptions:
- NAGIM samples still have dodgy headers (regex pattern avoids)
- Genomes, not Exomes
- 'to be avoided' sample list from other configs (may be out of date)
- pull SGs into analysis as whole family groups
"""

import re
from collections import defaultdict
from sys import argv

import toml

from metamist.graphql import gql, query


# CPG IDs matching this RE are probably NAGIM, so they need to be avoided
pat = re.compile(r'CPG2\d{5}')

# SG IDs marked as 'to be avoided' elsewhere in pipeline configs
# some of these may be out of date
REMOVE = [
    'CPG291633',
    'CPG291625',
    'CPG291617',
    'CPG291591',
    'CPG291609',
    'CPG291583',
    'CPG291567',
    'CPG291559',
    'CPG291575',
    'CPG11783',  # acute-care, no FASTQ data
    # 'CPG253328',  # perth-neuro, contamination rate 32% # put back in see https://github.com/populationgenomics/seqr-private/issues/41#issuecomment-1336537385
    'CPG13409',  # perth-neuro, coverage ~0x
    'CPG243717',
    # validation, NA12878_KCCG low coverage https://main-web.populationgenomics.org.au/validation/qc/cram/multiqc.html,
    'CPG246645',  # ag-hidden, eof issue  https://batch.hail.populationgenomics.org.au/batches/97645/jobs/440
    'CPG246678',  # ag-hidden, diff fastq size  https://batch.hail.populationgenomics.org.au/batches/97645/jobs/446
    'CPG246561',  # ag-hidden, coverage ~0x https://main-web.populationgenomics.org.au/ag-hidden/qc/cram/multiqc.html
    # rdp-kidney sequencing groups incorrectly imported into ibmdx
    'CPG261552',
    'CPG261560',
    'CPG261578',
    'CPG261586',
    'CPG261594',
    'CPG261602',
    'CPG261610',
    'CPG261628',
    'CPG261636',
    'CPG261644',
    'CPG261651',
    'CPG261669',
    'CPG261677',
    'CPG261685',
    'CPG261693',
    'CPG261701',
    'CPG261719',
    'CPG261727',
    'CPG261735',
    'CPG261743',
    'CPG261792',  # rdp-kidney misformated fastq - https://batch.hail.populationgenomics.org.au/batches/378736/jobs/43
    # acute care fasq parsing errors https://batch.hail.populationgenomics.org.au/batches/379303/jobs/24
    'CPG259150',
    'CPG258814',
    'CPG258137',
    'CPG258111',
    'CPG258012',
    # ohmr4 - sequencing groups that need to be deleted https://github.com/populationgenomics/seqr-private/issues/25
    'CPG221978',
    'CPG222000',
    # ohmr4 cram parsing in align issues
    'CPG261339',
    'CPG261347',
]

# skip these datasets, administrative/not RD
DULL_DATASETS = [
    'fewgenomes',
    'thousand-genomes',
    'tob-wgs',
    'hgdp',
    'seqr',
    'ibmdx',
    'mcri-lrp',
    'lof-curation',
    'talos',
    'udn-aus-training',
    'mito-mdt',
]

SEQUENCING_TYPE = 'genome'

SG_QUERY_STRING = """query SG_Query($project: String!){
  project(name: $project) {
    sequencingGroups {
      id
      type
      sample {
        participant {
          families {
            externalId
          }
        }
      }
    }
  }
}"""


def find_prior_analyses() -> set:
    """
    queries for all completed analyses, returns a set of SG IDs
    these SG IDs should not be picked up on subsequent runs

    Returns:
        set of all collected SGs
    """

    sv_analysis_query = gql(
        """
        query MyQuery {
            myProjects {
                analyses(type: {eq: "sv"}) {
                    meta
                    output
                }
            }
        }
    """
    )
    res = query(sv_analysis_query)
    sgs = set()
    for p in res['myProjects']:
        for a in p['analyses']:
            this_meta = a['meta']
            if (
                this_meta['stage'] != 'AnnotateVcf'
                or this_meta['sequencing_type'] != SEQUENCING_TYPE
            ):
                continue

            assert this_meta['type'] == 'gatk-sv-batch-calls', this_meta

            sgs.update(this_meta['sequencing_groups'])

    return sgs


def get_cpgs(project: str, dodge: set[str]) -> defaultdict:
    """
    Get the cpgs from the json, group on family Ext. ID
    """

    sgs_by_family = defaultdict(set)

    # query for results in this project
    this_json = query(SG_QUERY_STRING, variables={'project': project})

    for sg in this_json['project']['sequencingGroups']:
        if sg['type'] != SEQUENCING_TYPE:
            continue

        sgid = sg['id']

        # already analysed, or problematic
        if sgid in dodge or sgid in REMOVE:
            continue

        # limit to recent samples by ID, no NAGIM problems
        if not pat.match(sgid):
            continue

        # skip if no family
        try:
            ext = sg['sample']['participant']['families'][0]['externalId']
        except IndexError:
            continue
        sgs_by_family[ext].add(sgid)

    return sgs_by_family


def main(target_count: int, out_path: str):
    """
    Main function, queries for all projects, then gets the cpgs for each project
    Acquire SGs in family groups, until the min sample count is met

    Args:
        target_count (int): number of SGs we're aiming for
        out_path (str): write output list as JSON
    """

    # # query for all projects I can see
    # project_query = gql("""
    # query MyQuery {
    #     myProjects {
    #         dataset
    #     }
    # }
    # """)
    # for proj_data in query(project_query)['myProjects']:
    #     dataset_name = proj_data['dataset']

    # prepare something to catch all the SG IDs
    all_sgs_to_run: set[str] = set()
    datasets: set[str] = set()

    analysed_samples = find_prior_analyses()
    for dataset_name in [
        'udn-aus',
        'validation',
        'schr-neuro',
        'mito-disease',
        'ravenscroft-rdstudy',
        'ravenscroft-arch',
        'broad-rgp',
        'acute-care',
    ]:

        print(f'Processing {dataset_name}')

        # some duplicated, test projects?
        if dataset_name in DULL_DATASETS or dataset_name in datasets:
            print(f'Skipping {dataset_name}')
            continue

        datasets.add(dataset_name)

        # these are datasets we might be interested in
        dataset_sgs_by_family = get_cpgs(project=dataset_name, dodge=analysed_samples)

        print(f'Found {len(dataset_sgs_by_family)} families in {dataset_name}')

        # iterate until one or the other is exhausted
        while len(all_sgs_to_run) < target_count and dataset_sgs_by_family:
            # get the first family
            fam, sgs = dataset_sgs_by_family.popitem()

            # add the SGs to the total set
            all_sgs_to_run.update(sgs)

        # if the min count is not satisfied, pull next dataset
        if len(all_sgs_to_run) >= target_count:
            break

    with open(out_path, 'w', encoding='utf-8') as f:
        toml.dump(
            {
                'workflows': {
                    'input_datasets': list(datasets),
                    'only_sgs': list(all_sgs_to_run),
                }
            },
            f,
        )


if __name__ == '__main__':

    # maybe argparse/click this later
    # 2 args: min_samples, output file
    min_samples = int(argv[1])
    out_file = str(argv[2])
    main(min_samples, out_file)
