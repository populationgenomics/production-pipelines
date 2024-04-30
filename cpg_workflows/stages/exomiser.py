"""
As a dataset Stage

- Iterate over the families
- For each family export a VCF & subselected PED
- split the families into chunks
- for each chunk copy in the resources & set up a config file
- run exomiser, multiple families per setup
"""

from functools import cache

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_workflows.jobs.exomiser import (
    create_gvcf_to_vcf_jobs,
    extract_mini_ped_files,
    generate_seqr_summary,
    make_phenopackets,
    run_exomiser_13,
    run_exomiser_14,
)
from cpg_workflows.utils import get_logger
from cpg_workflows.workflow import Dataset, DatasetStage, SequencingGroup, StageInput, StageOutput, get_workflow, stage
from metamist.apis import ProjectApi

# this is used to separate family IDs from their individual outputs in RunExomiser
BREAKING_PUNCTUATION = '~~'
HPO_KEY: str = 'HPO Terms (present)'


@cache
def find_seqr_projects() -> dict[str, str]:
    """
    query for all the seqr projects
    map datasets to corresponding seqr names
    """

    project_api = ProjectApi()
    seq_type = config_retrieve(['workflow', 'sequencing_type'])
    seq_type_key = f'seqr-project-{seq_type}'
    return_dict: dict[str, str] = {}

    for element in project_api.get_seqr_projects():
        if element['meta'].get('is_seqr', False):
            if meta_proj := element['meta'].get(seq_type_key, ''):
                return_dict[element['dataset']] = meta_proj

    return return_dict


@cache
def find_families(dataset: Dataset) -> dict[str, list[SequencingGroup]]:
    """
    Find all the families in the project
    group on family ID and check for affected individuals & HPO terms
    re-group the selected members by an external ID
    at some point re-work this so that we have a better definition of proband
    """

    dict_by_family: dict[str, list[SequencingGroup]] = {}
    for sg in dataset.get_sequencing_groups():
        family_id = str(sg.pedigree.fam_id)
        dict_by_family.setdefault(family_id, []).append(sg)

    dict_by_ext_id: dict[str, list[SequencingGroup]] = {}
    # now remove families with no affected individuals
    for family, members in dict_by_family.items():

        # check for at least one retained member
        affected = [sg for sg in members if str(sg.pedigree.phenotype) == '2']

        # remove families with no affected members
        if not affected:
            get_logger(__file__).info(f'Family {family} has no affected individuals, skipping')
            continue

        # check that the affected members have HPO terms - required for exomiser
        if any([sg.meta['phenotypes'].get(HPO_KEY, '') == '' for sg in affected]):
            get_logger(__file__).info(f'Family {family} has affected individuals with no HPO terms, skipping')
            continue

        # key up those badbois using an affected external ID
        dict_by_ext_id[affected[0].external_id] = members

    return dict_by_ext_id


@stage
class CreateFamilyVCFs(DatasetStage):
    """
    it's a gVCF combiner, densification, Mt -> VCF
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        family_dict = find_families(dataset)
        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)
        return {
            str(family): dataset.prefix() / f'exomiser_{exomiser_version}_inputs' / f'{family}.vcf.bgz'
            for family in family_dict.keys()
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        family_dict = find_families(dataset)
        outputs = self.expected_outputs(dataset)
        jobs = create_gvcf_to_vcf_jobs(families=family_dict, out_paths=outputs)
        return self.make_outputs(dataset, outputs, jobs=jobs)


@stage
class MakePhenopackets(DatasetStage):
    """
    for each relevant family, make some Phenopackets
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)

        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)
        dataset_prefix = dataset.analysis_prefix() / f'exomiser_{exomiser_version}_inputs'
        return {family: dataset_prefix / f'{family}_phenopacket.json' for family in family_dict.keys()}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        this actually doesn't run as Jobs, but as a function...
        bit of an anti-pattern in this pipeline?
        """

        dataset_families = find_families(dataset)
        expected_out = self.expected_outputs(dataset)
        families_to_process = {k: v for k, v in dataset_families.items() if k in expected_out}
        make_phenopackets(families_to_process, expected_out)
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[])


@stage
class MakePedExtracts(DatasetStage):
    """
    from the dataset MT, we make a PED per-family
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)
        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)

        dataset_prefix = dataset.analysis_prefix() / f'exomiser_{exomiser_version}_inputs'
        return {family: dataset_prefix / f'{family}.ped' for family in family_dict.keys()}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        this actually doesn't run as Jobs, but as a function...
        bit of an anti-pattern in this pipeline?
        """
        dataset_families = find_families(dataset)
        expected_out = self.expected_outputs(dataset)
        families_to_process = {k: v for k, v in dataset_families.items() if k in expected_out}
        extract_mini_ped_files(families_to_process, expected_out)
        return self.make_outputs(dataset, data=expected_out, jobs=[])


@stage(required_stages=[CreateFamilyVCFs, MakePedExtracts, MakePhenopackets])
class RunExomiser(DatasetStage):
    """
    Run exomiser on the family VCFs

    Only run this for families with at least one affected individual
    """

    def expected_outputs(self, dataset: Dataset):
        """
        the pipeline logic only accepts a dictionary of paths
        but we need a dictionary of families, each containing a dictionary of paths
        I'm fudging this by putting some arbitrary punctuation in to split on later

        Args:
            dataset ():

        Returns:
            dict of outputs for this dataset
        """
        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)

        family_dict = find_families(dataset)
        dataset_prefix = dataset.analysis_prefix() / f'exomiser_{exomiser_version}_results'

        # only the TSVs are required
        return {family: dataset_prefix / f'{family}.tsv' for family in family_dict.keys()}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)

        output_dict = self.expected_outputs(dataset)

        vcf_inputs = inputs.as_dict(target=dataset, stage=CreateFamilyVCFs)
        ped_inputs = inputs.as_dict(target=dataset, stage=MakePedExtracts)
        ppk_files = inputs.as_dict(target=dataset, stage=MakePhenopackets)

        # combining all these stage outputs into one object
        single_dict = {
            family: {
                'output': output_dict[family],
                'vcf': vcf_inputs[family],
                'ped': ped_inputs[family],
                'pheno': ppk_files[family],
            }
            for family in output_dict.keys()
        }

        if exomiser_version == 14:
            jobs = run_exomiser_14(single_dict)
        elif exomiser_version == 13:
            jobs = run_exomiser_13(single_dict)
        else:
            raise ValueError(f'Exomiser version {exomiser_version} not supported')

        return self.make_outputs(dataset, data=output_dict, jobs=jobs)


@stage(required_stages=[RunExomiser], analysis_type='custom', analysis_keys=['tsv'])
class ExomiserSeqrTSV(DatasetStage):
    """
    Parse the Exomiser results into a TSV for Seqr
    """

    def expected_outputs(self, dataset: Dataset):
        exomiser_version = config_retrieve(['workflow', 'exomiser_version'], 14)

        return {
            'tsv': dataset.analysis_prefix()
            / get_workflow().output_version
            / f'exomiser_{exomiser_version}_results.tsv',
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        # is there a seqr project?
        projects = find_seqr_projects()
        if dataset.name not in projects:
            get_logger(__file__).info(f'No Seqr project found for {dataset.name}, skipping')
            return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[], skipped=True)

        results = inputs.as_dict(target=dataset, stage=RunExomiser)

        jobs = generate_seqr_summary(results, projects[dataset.name], str(self.expected_outputs(dataset)['tsv']))

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
