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
from cpg_workflows.jobs.exomiser import (
    create_gvcf_to_vcf_jobs,
    extract_mini_ped_files,
    make_phenopackets,
    run_exomiser_batches,
)
from cpg_workflows.utils import exists, get_logger
from cpg_workflows.workflow import Dataset, DatasetStage, SequencingGroup, StageInput, StageOutput, get_workflow, stage

BREAKING_PUNCTUATION = '~~'
HPO_KEY: str = 'HPO Terms (present)'


@cache
def find_families(dataset: Dataset) -> dict[str, list[SequencingGroup]]:
    """
    Find all the families in the project
    group on family ID and return the list of IDs
    """
    dict_by_family: dict[str, list[SequencingGroup]] = {}
    for sg in dataset.get_sequencing_groups():
        family_id = str(sg.pedigree.fam_id)
        dict_by_family.setdefault(family_id, []).append(sg)

    # now remove families with no affected individuals
    for family in list(dict_by_family.keys()):

        # check for at least one retained member
        affected = [sg for sg in dict_by_family[family] if str(sg.pedigree.phenotype) == '2']

        # remove families with no affected members
        if not affected:
            get_logger(__file__).info(f'Family {family} has no affected individuals, skipping')
            del dict_by_family[family]

        # # check that the affected members have HPO terms - required for exomiser
        # if any([sg.meta['phenotypes'].get(HPO_KEY, '') == '' for sg in affected]):
        #     get_logger(__file__).info(f'Family {family} has affected individuals with no HPO terms, skipping')
        #     del dict_by_family[family]
        #     continue

    return dict_by_family


@stage
class CreateFamilyVCFs(DatasetStage):
    """
    it's a gVCF combiner, densification, Mt -> VCF
    """

    def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
        family_dict = find_families(dataset)
        return {
            str(family): dataset.prefix() / 'exomiser_inputs' / f'{family}.vcf.bgz' for family in family_dict.keys()
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
        dataset_prefix = dataset.analysis_prefix() / 'exomiser_inputs'
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

        dataset_prefix = dataset.analysis_prefix() / 'exomiser_inputs'
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
        family_dict = find_families(dataset)
        dataset_prefix = dataset.analysis_prefix() / 'exomiser_results'

        # get all families without results
        output_files: dict[str, Path] = {}
        for family in family_dict.keys():

            # arbitrarily test one file
            if not exists(dataset_prefix / f'{family}_results.json'):
                output_files[f'{family}{BREAKING_PUNCTUATION}json'] = dataset_prefix / f'{family}_results.json'
                output_files[f'{family}{BREAKING_PUNCTUATION}gtsv'] = dataset_prefix / f'{family}_results.variants.tsv'
                output_files[f'{family}{BREAKING_PUNCTUATION}vtsv'] = dataset_prefix / f'{family}_results.genes.tsv'

        return output_files

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        output_dict = self.expected_outputs(dataset)

        # break out that dict again
        family_outputs: dict[str, dict[str, Path]] = {}
        for key, value in output_dict.items():
            family, file_type = key.split(BREAKING_PUNCTUATION)
            family_outputs.setdefault(family, {})[file_type] = value

        print(family_outputs)

        vcf_inputs = inputs.as_dict(target=dataset, stage=CreateFamilyVCFs)
        ped_inputs = inputs.as_dict(target=dataset, stage=MakePedExtracts)
        ppk_files = inputs.as_dict(target=dataset, stage=MakePhenopackets)

        # combining all these stage outputs into one object
        single_dict = {
            family: {
                'output': family_outputs[family],
                'vcf': vcf_inputs[family],
                'ped': ped_inputs[family],
                'pheno': ppk_files[family],
            }
            for family in family_outputs.keys()
        }

        jobs = run_exomiser_batches(single_dict)

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
