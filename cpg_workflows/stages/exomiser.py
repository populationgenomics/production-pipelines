"""
As a dataset Stage

- Iterate over the families
- For each family export a VCF & subselected PED
- split the families into chunks
- for each chunk copy in the resources & set up a config file
- run exomiser, multiple families per setup
"""

from functools import cache

from cpg_workflows.jobs.exomiser import (
    extract_mini_ped_files,
    extract_vcf_jobs,
    make_phenopackets,
    run_exomiser_batches,
)
from cpg_workflows.stages.aip import query_for_latest_mt
from cpg_workflows.utils import exists, get_logger
from cpg_workflows.workflow import Dataset, DatasetStage, SequencingGroup, StageInput, StageOutput, get_workflow, stage

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

        # check that the affected members have HPO terms - required for exomiser
        if any([sg.meta['phenotypes'].get(HPO_KEY, '') == '' for sg in affected]):
            get_logger(__file__).info(f'Family {family} has affected individuals with no HPO terms, skipping')
            del dict_by_family[family]
            continue

    return dict_by_family


# @stage
# class RDCombiner(DatasetStage):
#     """
#     it's a gVCF combiner, to be followed by densification
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
#         return {'vds': self.prefix / 'vds' / f'{get_workflow().output_version}.vds'}
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
#
#         sgids = dataset.get_sequencing_groups()
#         output = self.expected_outputs(dataset)
#         jobs = create_vds_jobs(sgids=sgids, out_path=str(output['vds']))
#         return self.make_outputs(dataset, output, jobs=jobs)
#
#
# @stage(required_stages=RDCombiner, analysis_keys=['mt'], analysis_type='custom')
# class VDStoMT(DatasetStage):
#     """
#     make a VDS from that MT
#     """
#
#     def expected_outputs(self, dataset: Dataset) -> dict[str, Path]:
#         return {'mt': self.prefix / 'mt' / f'{get_workflow().output_version}.mt'}
#
#     def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
#
#         vds = str(inputs.as_dict(target=dataset, stage=RDCombiner)['vds'])
#         output = self.expected_outputs(dataset)
#         jobs = mt_from_vds(vds_path=vds, out_path=str(output['mt']))
#         return self.make_outputs(dataset, output, jobs=jobs)


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
        # uses output_prefix to track consistently

        # the output path of the next step
        batch_prefix = get_workflow().prefix / 'RunExomiser'

        # get all families without results
        for family in family_dict.keys():
            if not exists(batch_prefix / f'{family}_results.json'):
                families_without_results.append(family)

        return {family: dataset.tmp_prefix() / self.name / f'{family}.vcf.bgz' for family in families_without_results}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:
        """ """
        dataset_families = find_families(dataset)
        expected_out = self.expected_outputs(dataset)
        families_to_process = {k: v for k, v in dataset_families.items() if k in expected_out}

        mt_path = query_for_latest_mt(dataset.name)
        vcf_jobs = extract_vcf_jobs(families_to_process, mt_path, expected_out)
        return self.make_outputs(dataset, data=expected_out, jobs=vcf_jobs)


@stage
class MakePhenopackets(DatasetStage):
    """
    for each relevant family, make some Phenopackets
    """

    def expected_outputs(self, dataset: Dataset):
        family_dict = find_families(dataset)
        families_without_results = []

        # the output path of the next step
        batch_prefix = get_workflow().prefix / 'RunExomiser'

        # get all families without results
        for family in family_dict.keys():
            if not exists(batch_prefix / f'{family}_results.json'):
                families_without_results.append(family)

        return {family: self.prefix / f'{family}_phenopacket.json' for family in families_without_results}

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
        families_without_results = []

        # get all families without results
        for family in family_dict.keys():
            if not exists(self.prefix / f'{family}_results.json'):
                families_without_results.append(family)

        return {family: self.prefix / f'{family}.ped' for family in families_without_results}

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
        family_dict = find_families(dataset)

        # only group samples with no result already
        # expecting this to fail sometimes
        families_without_results = []

        # get all families without results
        for family in family_dict.keys():
            if not exists(self.prefix / f'{family}_results.json'):
                families_without_results.append(family)
        return {family: self.prefix / f'{family}_results.json' for family in families_without_results}

    def queue_jobs(self, dataset: Dataset, inputs: StageInput) -> StageOutput:

        # todo filter out families with no affected individuals

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

        jobs = run_exomiser_batches(single_dict)

        return self.make_outputs(dataset, data=self.expected_outputs(dataset), jobs=jobs)
