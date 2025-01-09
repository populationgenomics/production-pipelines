"""
Targets for workflow stages: SequencingGroup, Dataset, Cohort.
"""

import copy
import hashlib
import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional

import pandas as pd

from cpg_utils import Path, to_path
from cpg_utils.config import dataset_path, get_config, reference_path, web_url

from .filetypes import AlignmentInput, BamPath, CramPath, FastqPairs, GvcfPath
from .metamist import Assay


def hash_from_list_of_strings(string_list: list[str], hash_length: int = 10, suffix: str | None = None) -> str:
    """
    Create a hash from a list of strings
    Args:
        string_list ():
        hash_length (int): how many characters to use from the hash
        suffix (str): optional, clarify the type of value which was hashed

    Returns:

    """
    hash_portion = hashlib.sha256(' '.join(string_list).encode()).hexdigest()[:hash_length]
    full_hash = f'{hash_portion}_{len(string_list)}'

    if suffix:
        full_hash += f'_{suffix}'
    return full_hash


class Target:
    """
    Defines a target that a stage can act upon.
    """

    def __init__(self) -> None:
        # Whether to process even if outputs exist:
        self.forced: bool = False
        # If not set, exclude from the workflow:
        self.active: bool = True

    def get_sequencing_groups(self, only_active: bool = True) -> list['SequencingGroup']:
        """
        Get flat list of all sequencing groups corresponding to this target.
        """
        raise NotImplementedError

    def get_sequencing_group_ids(self, only_active: bool = True) -> list[str]:
        """
        Get flat list of all sequencing group IDs corresponding to this target.
        """
        return [s.id for s in self.get_sequencing_groups(only_active=only_active)]

    def alignment_inputs_hash(self) -> str:
        """
        Unique hash string of sample alignment inputs. Useful to decide
        whether the analysis on the target needs to be rerun.
        """
        s = ' '.join(
            sorted(' '.join(str(s.alignment_input)) for s in self.get_sequencing_groups() if s.alignment_input),
        )
        h = hashlib.sha256(s.encode()).hexdigest()[:38]
        return f'{h}_{len(self.get_sequencing_group_ids())}'

    @property
    def target_id(self) -> str:
        """
        ID should be unique across target of all levels.

        We are raising NotImplementedError instead of making it an abstract class,
        because mypy is not happy about binding TypeVar to abstract classes, see:
        https://stackoverflow.com/questions/48349054/how-do-you-annotate-the-type-of
        -an-abstract-class-with-mypy

        Specifically,
        ```
        TypeVar('TargetT', bound=Target)
        ```
        Will raise:
        ```
        Only concrete class can be given where "Type[Target]" is expected
        ```
        """
        raise NotImplementedError

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        raise NotImplementedError

    def get_job_prefix(self) -> str:
        """
        Prefix job names.
        """
        raise NotImplementedError

    def rich_id_map(self) -> dict[str, str]:
        """
        Map if internal IDs to participant or external IDs, if the latter is provided.
        """
        return {s.id: s.rich_id for s in self.get_sequencing_groups() if s.participant_id != s.id}


class MultiCohort(Target):
    """
    Represents a "multi-cohort" target - multiple cohorts in the workflow.
    """

    def __init__(self) -> None:
        super().__init__()

        # previously MultiCohort.name was an underscore-delimited string of all the input cohorts
        # this was expanding to the point where filenames including this String were too long for *nix
        # instead we can create a hash of the input cohorts, and use that as the name
        # the exact cohorts can be obtained from the config associated with the ar-guid
        input_cohorts = get_config()['workflow'].get('input_cohorts', [])
        if input_cohorts:
            self.name = hash_from_list_of_strings(sorted(input_cohorts), suffix='cohorts')
        else:
            self.name = get_config()['workflow']['dataset']

        assert self.name, 'Ensure cohorts or dataset is defined in the config file.'

        self._cohorts_by_id: dict[str, Cohort] = {}
        self._datasets_by_name: dict[str, Dataset] = {}
        self.analysis_dataset = Dataset(name=get_config()['workflow']['dataset'])

    def __repr__(self):
        return f'MultiCohort({len(self.get_cohorts())} cohorts)'

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return self.name

    def get_cohorts(self, only_active: bool = True) -> list['Cohort']:
        """
        Gets list of all cohorts.
        Include only "active" cohorts (unless only_active is False)
        """
        cohorts = list(self._cohorts_by_id.values())
        if only_active:
            cohorts = [c for c in cohorts if c.active and c.get_datasets()]
        return cohorts

    def get_cohort_by_id(self, id: str, only_active: bool = True) -> Optional['Cohort']:
        """
        Get cohort by id.
        Include only "active" cohorts (unless only_active is False)
        """
        cohort = self._cohorts_by_id.get(id)
        if not cohort:
            logging.warning(f'Cohort {id} not found in the multi-cohort')
            return None
        if not only_active:  # Return cohort even if it's inactive
            return cohort
        if cohort.active and cohort.get_datasets():
            return cohort
        return None

    def get_datasets(self, only_active: bool = True) -> list['Dataset']:
        """
        Gets list of all datasets.
        Include only "active" datasets (unless only_active is False)
        """
        all_datasets = list(self._datasets_by_name.values())
        if only_active:
            all_datasets = [d for d in all_datasets if d.active and d.get_sequencing_groups()]
        return all_datasets

    def get_sequencing_groups(self, only_active: bool = True) -> list['SequencingGroup']:
        """
        Gets a flat list of all sequencing groups from all datasets.
        uses a dictionary to avoid duplicates (we could have the same sequencing group in multiple cohorts)
        Include only "active" sequencing groups (unless only_active is False)
        """
        all_sequencing_groups: dict[str, SequencingGroup] = {}
        for dataset in self.get_datasets(only_active):
            for sg in dataset.get_sequencing_groups(only_active):
                all_sequencing_groups[sg.id] = sg
        return list(all_sequencing_groups.values())

    def create_cohort(self, id: str, name: str) -> 'Cohort':
        """
        Create a cohort and add it to the multi-cohort.
        """
        if id in self._cohorts_by_id:
            logging.debug(f'Cohort {id} already exists in the multi-cohort')
            return self._cohorts_by_id[id]

        c = Cohort(id=id, name=name, multicohort=self)
        self._cohorts_by_id[c.id] = c
        return c

    def add_dataset(self, d: 'Dataset') -> 'Dataset':
        """
        Add a Dataset to the MultiCohort
        Args:
            d: Dataset object
        """
        if d.name in self._datasets_by_name:
            logging.debug(f'Dataset {d.name} already exists in the MultiCohort {self.name}')
        else:
            # We need create a new dataset to avoid manipulating the cohort dataset at this point
            self._datasets_by_name[d.name] = Dataset(d.name, d.cohort)
        return self._datasets_by_name[d.name]

    def get_dataset_by_name(self, name: str, only_active: bool = True) -> Optional['Dataset']:
        """
        Get dataset by name.
        Include only "active" datasets (unless only_active is False)
        """
        ds_by_name = {d.name: d for d in self.get_datasets(only_active)}
        return ds_by_name.get(name)

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        return {
            # 'sequencing_groups': self.get_sequencing_group_ids(),
            'datasets': [d.name for d in self.get_datasets()],
            'cohorts': [c.id for c in self.get_cohorts()],
        }

    def write_ped_file(self, out_path: Path | None = None, use_participant_id: bool = False) -> Path:
        """
        Create a PED file for all samples in the whole MultiCohort
        Duplication of the Cohort method
        PED is written with no header line to be strict specification compliant
        """
        datas = []
        for sequencing_group in self.get_sequencing_groups():
            datas.append(sequencing_group.pedigree.get_ped_dict(use_participant_id=use_participant_id))
        if not datas:
            raise ValueError(f'No pedigree data found for {self.name}')
        df = pd.DataFrame(datas)

        if out_path is None:
            out_path = self.analysis_dataset.tmp_prefix() / 'ped' / f'{self.alignment_inputs_hash()}.ped'

        if not get_config()['workflow'].get('dry_run', False):
            with out_path.open('w') as fp:
                df.to_csv(fp, sep='\t', index=False, header=False)
        return out_path


class Cohort(Target):
    """
    Represents a "cohort" target - all sequencing groups from a single CustomCohort (potentially spanning multiple datasets) in the workflow.
    Analysis dataset name is required and will be used as the default name for the
    cohort.
    """

    def __init__(self, id: str | None = None, name: str | None = None, multicohort: MultiCohort | None = None) -> None:
        super().__init__()
        self.id = id or get_config()['workflow']['dataset']
        self.name = name or get_config()['workflow']['dataset']
        self.analysis_dataset = Dataset(name=get_config()['workflow']['dataset'], cohort=self)
        self._datasets_by_name: dict[str, Dataset] = {}
        self.multicohort = multicohort

    def __repr__(self):
        return f'Cohort("{self.id}", {len(self.get_datasets())} datasets)'

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return self.id

    def write_ped_file(self, out_path: Path | None = None, use_participant_id: bool = False) -> Path:
        """
        Create a PED file for all samples in the whole cohort
        PED is written with no header line to be strict specification compliant
        """
        datas = []
        for sequencing_group in self.get_sequencing_groups():
            datas.append(sequencing_group.pedigree.get_ped_dict(use_participant_id=use_participant_id))
        if not datas:
            raise ValueError(f'No pedigree data found for {self.id}')
        df = pd.DataFrame(datas)

        if out_path is None:
            out_path = self.analysis_dataset.tmp_prefix() / 'ped' / f'{self.alignment_inputs_hash()}.ped'

        if not get_config()['workflow'].get('dry_run', False):
            with out_path.open('w') as fp:
                df.to_csv(fp, sep='\t', index=False, header=False)
        return out_path

    def get_datasets(self, only_active: bool = True) -> list['Dataset']:
        """
        Gets list of all datasets.
        Include only "active" datasets (unless only_active is False)
        """
        datasets = [ds for k, ds in self._datasets_by_name.items()]
        if only_active:
            datasets = [ds for ds in datasets if ds.active and ds.get_sequencing_groups()]
        return datasets

    def get_dataset_by_name(self, name: str, only_active: bool = True) -> Optional['Dataset']:
        """
        Get dataset by name.
        Include only "active" datasets (unless only_active is False)
        """
        ds_by_name = {d.name: d for d in self.get_datasets(only_active)}
        return ds_by_name.get(name)

    def get_sequencing_groups(self, only_active: bool = True) -> list['SequencingGroup']:
        """
        Gets a flat list of all sequencing groups from all datasets.
        Include only "active" sequencing groups (unless only_active is False)
        """
        all_sequencing_groups = []
        for ds in self.get_datasets(only_active=False):
            all_sequencing_groups.extend(ds.get_sequencing_groups(only_active=only_active))
        return all_sequencing_groups

    def add_dataset(self, dataset: 'Dataset') -> 'Dataset':
        """
        Add existing dataset into the cohort.
        """
        dataset.cohort = self
        if dataset.name in self._datasets_by_name:
            logging.debug(f'Dataset {dataset.name} already exists in the cohort')
            return dataset
        self._datasets_by_name[dataset.name] = dataset
        return dataset

    def create_dataset(self, name: str) -> 'Dataset':
        """
        Create a dataset and add it to the cohort.
        """
        if name in self._datasets_by_name:
            logging.debug(f'Dataset {name} already exists in the cohort')
            return self._datasets_by_name[name]

        if name == self.analysis_dataset.name:
            ds = self.analysis_dataset
        else:
            ds = Dataset(name=name, cohort=self)

        self._datasets_by_name[ds.name] = ds
        return ds

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        return {
            # 'sequencing_groups': self.get_sequencing_group_ids(),
            'datasets': [d.name for d in self.get_datasets()],
        }

    def get_job_prefix(self) -> str:
        """
        Prefix job names.
        """
        return ''

    def to_tsv(self) -> str:
        """
        Export to a parsable TSV file
        """
        assert self.get_sequencing_groups()
        tsv_path = self.analysis_dataset.tmp_prefix() / 'samples.tsv'
        df = pd.DataFrame(
            {
                's': s.id,
                'gvcf': s.gvcf or '-',
                'sex': s.meta.get('sex') or '-',
                'continental_pop': s.meta.get('continental_pop') or '-',
                'subcontinental_pop': s.meta.get('subcontinental_pop') or '-',
            }
            for s in self.get_sequencing_groups()
        ).set_index('s', drop=False)
        with to_path(tsv_path).open('w') as f:
            df.to_csv(f, index=False, sep='\t', na_rep='NA')
        return tsv_path


class Dataset(Target):
    """
    Represents a CPG dataset.

    Each `dataset` at the CPG corresponds to
    * a GCP project: https://github.com/populationgenomics/team-docs/tree/main/storage_policies
    * a Pulumi stack: https://github.com/populationgenomics/analysis-runner/tree/main/stack
    * a metamist project
    """

    def __init__(
        self,
        name: str,
        cohort: Cohort | None = None,
    ):
        super().__init__()
        self._sequencing_group_by_id: dict[str, SequencingGroup] = {}
        self.name = name
        self.cohort = cohort
        self.active = True

    @staticmethod
    def create(name: str) -> 'Dataset':
        """
        Create a dataset.
        """
        return Dataset(name=name)

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return self.name

    def __repr__(self):
        return f'Dataset("{self.name}", {len(self.get_sequencing_groups())} sequencing groups)'

    def __str__(self):
        return f'{self.name} ({len(self.get_sequencing_groups())} sequencing groups)'

    def _seq_type_subdir(self) -> str:
        """
        Subdirectory parametrised by sequencing type. For genomes, we don't
        prefix at all.
        """
        seq_type = get_config()['workflow'].get('sequencing_type')
        return '' if not self.cohort or not seq_type or seq_type == 'genome' else seq_type

    def prefix(self, **kwargs) -> Path:
        """
        The primary storage path.
        """
        return to_path(
            dataset_path(
                self._seq_type_subdir(),
                dataset=self.name,
                **kwargs,
            ),
        )

    def tmp_prefix(self, **kwargs) -> Path:
        """
        Storage path for temporary files.
        """
        return to_path(
            dataset_path(
                self._seq_type_subdir(),
                dataset=self.name,
                category='tmp',
                **kwargs,
            ),
        )

    def analysis_prefix(self, **kwargs) -> Path:
        """
        Storage path for analysis files.
        """
        return to_path(
            dataset_path(
                self._seq_type_subdir(),
                dataset=self.name,
                category='analysis',
                **kwargs,
            ),
        )

    def web_prefix(self, **kwargs) -> Path:
        """
        Path for files served by an HTTP server Matches corresponding URLs returns by
        self.web_url() URLs.
        """
        return to_path(
            dataset_path(
                self._seq_type_subdir(),
                dataset=self.name,
                category='web',
                **kwargs,
            ),
        )

    def web_url(self) -> str | None:
        """
        URLs matching self.storage_web_path() files serverd by an HTTP server.
        """
        return web_url(
            self._seq_type_subdir(),
            dataset=self.name,
        )

    def add_sequencing_group(
        self,
        id: str,  # pylint: disable=redefined-builtin
        *,
        sequencing_type: str,
        sequencing_technology: str,
        sequencing_platform: str,
        external_id: str | None = None,
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Optional['Sex'] = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input: AlignmentInput | None = None,
    ) -> 'SequencingGroup':
        """
        Create a new sequencing group and add it to the dataset.
        """
        if id in self._sequencing_group_by_id:
            logging.debug(f'SequencingGroup {id} already exists in the dataset {self.name}')
            return self._sequencing_group_by_id[id]

        force_sgs = get_config()['workflow'].get('force_sgs', set())
        forced = id in force_sgs or external_id in force_sgs or participant_id in force_sgs

        s = SequencingGroup(
            id=id,
            dataset=self,
            external_id=external_id,
            sequencing_type=sequencing_type,
            sequencing_technology=sequencing_technology,
            sequencing_platform=sequencing_platform,
            participant_id=participant_id,
            meta=meta,
            sex=sex,
            pedigree=pedigree,
            alignment_input=alignment_input,
            forced=forced,
        )
        self._sequencing_group_by_id[id] = s
        return s

    def add_sequencing_group_object(self, s: 'SequencingGroup'):
        """
        Add a sequencing group object to the dataset.
        Args:
            s: SequencingGroup object
        """
        if s.id in self._sequencing_group_by_id:
            logging.debug(f'SequencingGroup {s.id} already exists in the dataset {self.name}')
            return self._sequencing_group_by_id[s.id]
        self._sequencing_group_by_id[s.id] = s

    def get_sequencing_groups(self, only_active: bool = True) -> list['SequencingGroup']:
        """
        Get dataset's sequencing groups. Include only "active" sequencing groups, unless only_active=False
        """
        return [s for sid, s in self._sequencing_group_by_id.items() if (s.active or not only_active)]

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        return {
            'dataset': self.name,
            # 'sequencing_groups': self.get_sequencing_group_ids(),
        }

    def get_job_prefix(self) -> str:
        """
        Prefix job names.
        """
        return f'{self.name}: '

    def write_ped_file(self, out_path: Path | None = None, use_participant_id: bool = False) -> Path:
        """
        Create a PED file for all sequencing groups
        PED is written with no header line to be strict specification compliant
        """
        datas = []
        for sequencing_group in self.get_sequencing_groups():
            datas.append(sequencing_group.pedigree.get_ped_dict(use_participant_id=use_participant_id))
        if not datas:
            raise ValueError(f'No pedigree data found for {self.name}')
        df = pd.DataFrame(datas)

        if out_path is None:
            out_path = self.tmp_prefix() / 'ped' / f'{self.alignment_inputs_hash()}.ped'

        if not get_config()['workflow'].get('dry_run', False):
            with out_path.open('w') as fp:
                df.to_csv(fp, sep='\t', index=False, header=False)
        return out_path


class Sex(Enum):
    """
    Sex as in PED format
    """

    UNKNOWN = 0
    MALE = 1
    FEMALE = 2

    @staticmethod
    def parse(sex: str | int | None) -> 'Sex':
        """
        Parse a string into a Sex object.
        """
        if sex:
            _sex = sex.lower() if isinstance(sex, str) else sex
            if _sex in ('m', 'male', '1', 1):
                return Sex.MALE
            if _sex in ('f', 'female', '2', 2):
                return Sex.FEMALE
            if _sex in ('u', 'unknown', '0', 0):
                return Sex.UNKNOWN
            raise ValueError(f'Unrecognised sex value {sex}')
        return Sex.UNKNOWN

    def __str__(self):
        return self.name


class SequencingGroup(Target):
    """
    Represents a sequencing group.
    """

    def __init__(
        self,
        id: str,  # pylint: disable=redefined-builtin
        dataset: 'Dataset',  # type: ignore  # noqa: F821
        *,
        sequencing_type: str,
        sequencing_technology: str,
        sequencing_platform: str,
        external_id: str | None = None,
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Sex | None = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input: AlignmentInput | None = None,
        assays: tuple[Assay, ...] | None = None,
        forced: bool = False,
    ):
        super().__init__()
        self.id = id
        self.name = id
        self._external_id = external_id
        self.sequencing_type = sequencing_type
        self.sequencing_technology = sequencing_technology
        self.sequencing_platform = sequencing_platform

        self.dataset = dataset
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.pedigree: PedigreeInfo = pedigree or PedigreeInfo(
            sequencing_group=self,
            fam_id=self.participant_id,
            sex=sex or Sex.UNKNOWN,
        )
        if sex:
            self.pedigree.sex = sex
        self.alignment_input: AlignmentInput | None = alignment_input
        self.assays: tuple[Assay, ...] | None = assays
        self.forced = forced
        self.active = True
        # Only set if the file exists / found in Metamist:
        self.gvcf: GvcfPath | None = None
        self.cram: CramPath | None = None

    def __repr__(self):
        values = {
            'participant': self._participant_id if self._participant_id else '',
            'sequencing_type': self.sequencing_type,
            'sequencing_technology': self.sequencing_technology,
            'sequencing_platform': self.sequencing_platform,
            'forced': str(self.forced),
            'active': str(self.active),
            'meta': str(self.meta),
            'alignment_inputs': self.alignment_input,
            'pedigree': self.pedigree,
        }
        retval = f'SequencingGroup({self.dataset.name}/{self.id}'
        if self._external_id:
            retval += f'|{self._external_id}'
        return retval + ''.join(f', {k}={v}' for k, v in values.items())

    def __str__(self):
        ai_tag = ''
        if self.alignment_input:
            ai_tag += f'|SEQ={self.sequencing_type}:'
            if isinstance(self.alignment_input, CramPath):
                ai_tag += 'CRAM'
            elif isinstance(self.alignment_input, BamPath):
                ai_tag += 'BAM'
            elif isinstance(self.alignment_input, FastqPairs):
                ai_tag += f'{len(self.alignment_input)}FQS'
            else:
                raise ValueError(f'Unrecognised alignment input type {type(self.alignment_input)}')

        ext_id = f'|{self._external_id}' if self._external_id else ''
        return f'SequencingGroup({self.dataset.name}/{self.id}{ext_id}{ai_tag})'

    @property
    def participant_id(self) -> str:
        """
        Get ID of participant corresponding to this sequencing group,
        or substitute it with external ID.
        """
        return self._participant_id or self.external_id

    @participant_id.setter
    def participant_id(self, val: str):
        """
        Set participant ID.
        """
        self._participant_id = val

    @property
    def external_id(self) -> str:
        """
        Get external sample ID, or substitute it with the internal ID.
        """
        return self._external_id or self.id

    @property
    def rich_id(self) -> str:
        """
        ID for reporting purposes: composed of internal as well as external
        or participant IDs.
        """
        return self.id + '|' + self.participant_id

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sequencing group, corresponding
        a PED file entry.
        """
        return self.pedigree.get_ped_dict(use_participant_id)

    def make_cram_path(self) -> CramPath:
        """
        Path to a CRAM file. Not checking its existence here.
        """
        path = self.dataset.prefix() / 'cram' / f'{self.id}.cram'
        return CramPath(
            path=path,
            index_path=path.with_suffix('.cram.crai'),
            reference_assembly=reference_path('broad/ref_fasta'),
        )

    def make_gvcf_path(self) -> GvcfPath:
        """
        Path to a GVCF file. Not checking its existence here.
        """
        return GvcfPath(self.dataset.prefix() / 'gvcf' / f'{self.id}.g.vcf.gz')

    @property
    def make_sv_evidence_path(self) -> Path:
        """
        Path to the evidence root for GATK-SV evidence files.
        """
        return self.dataset.prefix() / 'sv_evidence'

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return self.id

    def get_sequencing_groups(self, only_active: bool = True) -> list['SequencingGroup']:
        """
        Implementing the abstract method.
        """
        if only_active and not self.active:
            return []
        return [self]

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        attrs = {
            'dataset': self.dataset.name,
            'sequencing_group': self.id,
        }
        _participant_id: str | None = self._participant_id or self._external_id
        if _participant_id:
            attrs['participant_id'] = _participant_id
        return attrs

    def get_job_prefix(self) -> str:
        """
        Prefix job names.
        """
        return f'{self.dataset.name}/{self.id}: '


@dataclass
class PedigreeInfo:
    """
    Pedigree relationships with other sequencing groups in the cohort, and other PED data
    """

    sequencing_group: SequencingGroup
    sex: Sex = Sex.UNKNOWN
    fam_id: str | None = None
    phenotype: str | int = 0
    dad: SequencingGroup | None = None
    mom: SequencingGroup | None = None

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sequencing group, corresponding
        a PED file entry.
        """

        def _get_id(_s: SequencingGroup | None) -> str:
            if _s is None:
                return '0'
            if use_participant_id:
                return _s.participant_id
            return _s.id

        return {
            'Family.ID': self.fam_id or self.sequencing_group.participant_id,
            'Individual.ID': _get_id(self.sequencing_group),
            'Father.ID': _get_id(self.dad),
            'Mother.ID': _get_id(self.mom),
            'Sex': str(self.sex.value),
            'Phenotype': str(self.phenotype),
        }
