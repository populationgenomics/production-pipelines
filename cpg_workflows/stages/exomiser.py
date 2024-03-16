"""
As a dataset Stage

- Iterate over the families
- For each family export a VCF & subselected PED
- split the families into chunks
- for each chunk copy in the resources & set up a config file
- run exomiser, multiple families per setup
"""

from functools import lru_cache

from cpg_workflows.jobs.exomiser import (
    extract_vcf_jobs,
    extract_mini_ped_files,
    run_exomiser_batches,
)
from cpg_workflows.stages.aip import query_for_latest_mt
from cpg_workflows.workflow import (
    get_workflow,
    StageInput,
    StageOutput,
    DatasetStage,
    Dataset,
    stage,
)
from cpg_workflows.utils import exists


@lru_cache(maxsize=0)
def find_families(dataset: Dataset) -> dict[str, list[str]]:
    """
    Find all the families in the project
    group on family ID and return the list of IDs
    """
    dict_by_family: dict[str, list[str]] = {}
    for sg in dataset.get_sequencing_groups():
        family_id = sg.pedigree.fam_id
        assert isinstance(family_id, str), f'Family ID is not a string: {family_id}'
        dict_by_family.setdefault(family_id, []).append(sg.id)

    return dict_by_family


@stage
class CreateFamilyVCFs(DatasetStage):
    """
    from the dataset MT, we make a VCF per-family
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)
        families_without_results = []

        # generate the exomiser result paths for each family
        # skip those families that already have results
        # uses output_prefix to track consisently

        # the output path of the next step
        batch_prefix = get_workflow().prefix / RunExomiser.name

        # get all families without results
        for family in family_dict.keys():
            if not exists(batch_prefix / f'{family}_results.json'):
                families_without_results.append(family)

        return {
            family: dataset.tmp_prefix() / self.name / f'{family}.vcf.bgz'
            for family in families_without_results
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """ """
        dataset_families = find_families(dataset)
        expected_out = self.expected_outputs(dataset)
        families_to_process = {
            k: v for k, v in dataset_families.items() if k in expected_out
        }
        vcf_jobs = extract_vcf_jobs(
            families_to_process,
            query_for_latest_mt(dataset.name),
            dataset.tmp_prefix() / self.name,
        )
        return self.make_outputs(dataset, data=expected_out, jobs=vcf_jobs)


@stage
class MakePedExtracts(DatasetStage):
    """
    from the dataset MT, we make a PED per-family
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)
        families_without_results = []

        # get all families without results
        for family in family_dict.keys():
            if not exists(self.prefix / f'{family}_results.json'):
                families_without_results.append(family)

        return {
            family: self.prefix / f'{family}.ped' for family in families_without_results
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """
        this actually doesn't run as Jobs, but as a function...
        bit of an anti-pattern in this pipeline?
        """
        dataset_families = find_families(dataset)
        expected_out = self.expected_outputs(dataset)
        families_to_process = {
            k: v for k, v in dataset_families.items() if k in expected_out
        }
        out_path = dataset.tmp_prefix() / self.name
        extract_mini_ped_files(dataset, families_to_process, out_path)
        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[])


@stage(required_stages=[CreateFamilyVCFs, MakePedExtracts])
class RunExomiser(DatasetStage):
    """
    Run exomiser on the family VCFs, using the mini-vcf and mini-ped
    HUGELY anti-pattern'y, this is a DatasetStage, but it really just
    captures SequencingGroups in chunks to optimise the exomiser jobs
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)

        # only group samples with no result already
        # expecting this to fail sometimes
        families_without_results = []

        # get all families without results
        for family in family_dict.keys():
            if not exists(self.prefix / f'{family}_results.json'):
                families_without_results.append(family)
        return {
            family: dataset.prefix() / f'{family}_results.json'
            for family in family_dict.keys()
        }

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        output_dict = self.expected_outputs(dataset)
        vcf_inputs = inputs.as_dict(target=dataset, stage=CreateFamilyVCFs)
        ped_inputs = inputs.as_dict(target=dataset, stage=MakePedExtracts)

        # expecting all the keys to be the same...
        single_dict = {
            family: {
                'output': output_dict[family],
                'vcf': vcf_inputs[family],
                'ped': ped_inputs[family],
            }
            for family in output_dict.keys()
        }

        jobs = run_exomiser_batches(single_dict)

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=[])
