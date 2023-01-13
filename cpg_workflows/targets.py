"""
Targets for workflow stages: Sample, Dataset, Cohort.
"""

import hashlib
import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import pandas as pd

from cpg_utils.hail_batch import dataset_path, web_url, reference_path
from cpg_utils.config import get_config
from cpg_utils import Path, to_path

from .filetypes import (
    AlignmentInput,
    CramPath,
    BamPath,
    GvcfPath,
    FastqPairs,
)
from .metamist import Sequence


class Target:
    """
    Defines a target that a stage can act upon.
    """

    def __init__(self):
        # Whether to process even if outputs exist:
        self.forced: bool = False
        # If not set, exclude from the workflow:
        self.active: bool = True

    def get_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Get flat list of all samples corresponding to this target.
        """
        raise NotImplementedError

    def get_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Get flat list of all sample IDs corresponding to this target.
        """
        return [s.id for s in self.get_samples(only_active=only_active)]

    def alignment_inputs_hash(self) -> str:
        """
        Unique hash string of sample alignment inputs. Useful to decide
        whether the analysis on the target needs to be rerun.
        """
        s = ' '.join(
            sorted(
                [
                    ' '.join(
                        sorted(
                            str(alignment_input)
                            for alignment_input in s.alignment_input_by_seq_type.values()
                        )
                    )
                    for s in self.get_samples()
                    if s.alignment_input_by_seq_type
                ]
            )
        )
        h = hashlib.sha256(s.encode()).hexdigest()[:38]
        return f'{h}_{len(self.get_sample_ids())}'

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
        return {s.id: s.rich_id for s in self.get_samples() if s.participant_id != s.id}


class Cohort(Target):
    """
    Represents a "cohort" target - all samples from all datasets in the workflow.
    Analysis dataset name is required and will be used as the default name for the
    cohort.
    """

    def __init__(self):
        super().__init__()
        self.name = get_config()['workflow']['dataset']
        self.analysis_dataset = Dataset(name=self.name, cohort=self)
        self._datasets_by_name: dict[str, Dataset] = {}

    def __repr__(self):
        return f'Cohort("{self.name}", {len(self.get_datasets())} datasets)'

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return self.name

    def get_datasets(self, only_active: bool = True) -> list['Dataset']:
        """
        Gets list of all datasets.
        Include only "active" datasets (unless only_active is False)
        """
        datasets = [ds for k, ds in self._datasets_by_name.items()]
        if only_active:
            datasets = [ds for ds in datasets if ds.active and ds.get_samples()]
        return datasets

    def get_dataset_by_name(
        self, name: str, only_active: bool = True
    ) -> Optional['Dataset']:
        """
        Get dataset by name.
        Include only "active" datasets (unless only_active is False)
        """
        ds_by_name = {d.name: d for d in self.get_datasets(only_active)}
        return ds_by_name.get(name)

    def get_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Gets a flat list of all samples from all datasets.
        Include only "active" samples (unless only_active is False)
        """
        all_samples = []
        for ds in self.get_datasets(only_active=False):
            all_samples.extend(ds.get_samples(only_active=only_active))
        return all_samples

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

    def create_dataset(
        self,
        name: str,
    ) -> 'Dataset':
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
            'samples': self.get_sample_ids(),
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
        assert self.get_samples()
        tsv_path = self.analysis_dataset.tmp_prefix() / 'samples.tsv'
        df = pd.DataFrame(
            {
                's': s.id,
                'gvcf': s.gvcf or '-',
                'sex': s.meta.get('sex') or '-',
                'continental_pop': s.meta.get('continental_pop') or '-',
                'subcontinental_pop': s.meta.get('subcontinental_pop') or '-',
            }
            for s in self.get_samples()
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
        self._sample_by_id: dict[str, Sample] = {}
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
        return f'Dataset("{self.name}", {len(self.get_samples())} samples)'

    def __str__(self):
        return f'{self.name} ({len(self.get_samples())} samples)'

    def _seq_type_subdir(self) -> str:
        """
        Subdirectory parametrised by sequencing type. For genomes, we don't
        prefix at all.
        """
        seq_type = get_config()['workflow'].get('sequencing_type')
        return (
            '' if not self.cohort or not seq_type or seq_type == 'genome' else seq_type
        )

    def prefix(self, **kwargs) -> Path:
        """
        The primary storage path.
        """
        return to_path(
            dataset_path(
                self._seq_type_subdir(),
                dataset=self.name,
                **kwargs,
            )
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
            )
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
            )
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
            )
        )

    def web_url(self) -> str | None:
        """
        URLs matching self.storage_web_path() files serverd by an HTTP server.
        """
        return web_url(
            self._seq_type_subdir(),
            dataset=self.name,
        )

    def add_sample(
        self,
        id: str,  # pylint: disable=redefined-builtin
        external_id: str | None = None,
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Optional['Sex'] = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input_by_seq_type: dict[str, AlignmentInput] | None = None,
    ) -> 'Sample':
        """
        Create a new sample and add it to the dataset.
        """
        if id in self._sample_by_id:
            logging.debug(f'Sample {id} already exists in the dataset {self.name}')
            return self._sample_by_id[id]

        force_samples = get_config()['workflow'].get('force_samples', set())
        forced = (
            id in force_samples
            or external_id in force_samples
            or participant_id in force_samples
        )

        s = Sample(
            id=id,
            dataset=self,
            external_id=external_id,
            participant_id=participant_id,
            meta=meta,
            sex=sex,
            pedigree=pedigree,
            alignment_input_by_seq_type=alignment_input_by_seq_type,
            forced=forced,
        )
        self._sample_by_id[id] = s
        return s

    def get_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Get dataset's samples. Include only "active" samples, unless only_active=False
        """
        return [
            s for sid, s in self._sample_by_id.items() if (s.active or not only_active)
        ]

    def get_job_attrs(self) -> dict:
        """
        Attributes for Hail Batch job.
        """
        return {
            'dataset': self.name,
            'samples': self.get_sample_ids(),
        }

    def get_job_prefix(self) -> str:
        """
        Prefix job names.
        """
        return f'{self.name}: '

    def write_ped_file(
        self, out_path: Path | None = None, use_participant_id: bool = False
    ) -> Path:
        """
        Create a PED file for all samples
        """
        datas = []
        for sample in self.get_samples():
            datas.append(
                sample.pedigree.get_ped_dict(use_participant_id=use_participant_id)
            )
        if not datas:
            raise ValueError(f'No pedigree data found for {self.name}')
        df = pd.DataFrame(datas)

        if out_path is None:
            out_path = self.tmp_prefix() / 'ped' / f'{self.alignment_inputs_hash()}.ped'

        if not get_config()['workflow'].get('dry_run', False):
            with out_path.open('w') as fp:
                df.to_csv(fp, sep='\t', index=False)
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


class Sample(Target):
    """
    Represents a Sample.
    """

    def __init__(
        self,
        id: str,  # pylint: disable=redefined-builtin
        dataset: 'Dataset',  # type: ignore  # noqa: F821
        external_id: str | None = None,
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Sex | None = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input_by_seq_type: dict[str, AlignmentInput] | None = None,
        seq_by_type: dict[str, Sequence] | None = None,
        forced: bool = False,
    ):
        super().__init__()
        self.id = id
        self._external_id = external_id
        self.dataset = dataset
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.pedigree: PedigreeInfo = pedigree or PedigreeInfo(
            sample=self,
            fam_id=self.participant_id,
            sex=sex or Sex.UNKNOWN,
        )
        if sex:
            self.pedigree.sex = sex
        self.alignment_input_by_seq_type: dict[str, AlignmentInput] = (
            alignment_input_by_seq_type or dict()
        )
        self.seq_by_type: dict[str, Sequence] = seq_by_type or dict()
        self.forced = forced
        self.active = True
        # Only set if the file exists / found in Metamist:
        self.gvcf: GvcfPath | None = None
        self.cram: CramPath | None = None

    def __repr__(self):
        values = {
            'participant': self._participant_id if self._participant_id else '',
            'forced': str(self.forced),
            'active': str(self.active),
            'meta': str(self.meta),
            'alignment_inputs': ','.join(
                [
                    f'{seq_t}: {al_inp}'
                    for seq_t, al_inp in self.alignment_input_by_seq_type.items()
                ]
            ),
            'pedigree': self.pedigree,
        }
        retval = f'Sample({self.dataset.name}/{self.id}'
        if self._external_id:
            retval += f'|{self._external_id}'
        return retval + ''.join(f', {k}={v}' for k, v in values.items())

    def __str__(self):
        ai_tag = ''
        for seq_type, alignment_input in self.alignment_input_by_seq_type.items():
            ai_tag += f'|SEQ={seq_type}:'
            if isinstance(alignment_input, CramPath):
                ai_tag += 'CRAM'
            elif isinstance(alignment_input, BamPath):
                ai_tag += 'BAM'
            else:
                assert isinstance(alignment_input, FastqPairs)
                ai_tag += f'{len(alignment_input)}FQS'

        ext_id = f'|{self._external_id}' if self._external_id else ''
        return f'Sample({self.dataset.name}/{self.id}{ext_id}{ai_tag})'

    @property
    def participant_id(self) -> str:
        """
        Get ID of participant corresponding to this sample,
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
        Returns a dictionary of pedigree fields for this sample, corresponding
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
    def target_id(self) -> str:
        """Unique target ID"""
        return self.id

    def get_samples(self, only_active: bool = True) -> list['Sample']:
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
            'sample': self.id,
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
    Pedigree relationships with other samples in the cohort, and other PED data
    """

    sample: Sample
    sex: Sex = Sex.UNKNOWN
    fam_id: str | None = None
    phenotype: str | int = 0
    dad: Sample | None = None
    mom: Sample | None = None

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sample, corresponding
        a PED file entry.
        """

        def _get_id(_s: Sample | None) -> str:
            if _s is None:
                return '0'
            if use_participant_id:
                return _s.participant_id
            return _s.id

        return {
            'Family.ID': self.fam_id or self.sample.participant_id,
            'Individual.ID': _get_id(self.sample),
            'Father.ID': _get_id(self.dad),
            'Mother.ID': _get_id(self.mom),
            'Sex': str(self.sex.value),
            'Phenotype': str(self.phenotype),
        }
