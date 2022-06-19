"""
Pipeline stage target and subclasses: Cohort, Dataset, Sample.
Plus Sample properties: PedigreeInfo, Sex.
"""
import hashlib
import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional
import pandas as pd
from cpg_utils.hail_batch import dataset_path, web_url

from . import Namespace, Path, to_path
from .types import AlignmentInput, CramPath, GvcfPath, SequencingType, FastqPairs

logger = logging.getLogger(__file__)


class Target:
    """
    Defines a target that stage can act upon. Classes like Sample, Dataset, Pipeline
    extend this class.
    """

    def __init__(self):
        # Whether to process even if outputs exist:
        self.forced: bool = False
        # If not set, exclude from the pipeline:
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
                    ' '.join(sorted(
                        str(alignment_input)
                        for alignment_input in s.alignment_input_by_seq_type.values()
                    ))
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

    def external_id_map(self) -> dict[str, str]:
        """
        Map if internal IDs to participant or external IDs, if the latter is provided.
        """
        return {
            s.id: s.id + '|' + s.participant_id
            for s in self.get_samples()
            if s.participant_id != s.id
        }


class Cohort(Target):
    """
    Represents a "cohort" target - all samples from all datasets in the pipeline.
    Analysis dataset name is required and will be used as the default name for the
    cohort.
    """

    def __init__(
        self,
        analysis_dataset_name: str,
        namespace: Namespace,
        name: str | None = None,
        sequencing_type: SequencingType = SequencingType.GENOME,
    ):
        super().__init__()
        self.name = name or analysis_dataset_name
        self.namespace = namespace
        self.analysis_dataset = Dataset(
            name=analysis_dataset_name,
            namespace=namespace,
            cohort=self,
        )
        self.sequencing_type = sequencing_type
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
            logger.debug(f'Dataset {dataset.name} already exists in the cohort')
            return dataset
        self._datasets_by_name[dataset.name] = dataset
        return dataset

    def create_dataset(
        self,
        name: str,
        namespace: Namespace | None = None,
    ) -> 'Dataset':
        """
        Create a dataset and add it to the cohort.
        """
        namespace = namespace or self.analysis_dataset.namespace
        # Normalising the dataset's name:
        name = build_dataset_name(*parse_stack(name, namespace))
        if name in self._datasets_by_name:
            logger.debug(f'Dataset {name} already exists in the cohort')
            return self._datasets_by_name[name]

        if name == self.analysis_dataset.name:
            ds = self.analysis_dataset
        else:
            ds = Dataset(
                name=name,
                namespace=namespace,
                cohort=self,
            )

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


def parse_stack(name: str, namespace: Namespace | None = None) -> tuple[str, Namespace]:
    """
    Input `name` can be either e.g. "seqr" or "seqr-test". The latter will be
    resolved to stack="seqr" and is_test=True, unless `namespace` is provided
    explicitly.

    Returns the stack id and a corrected namespace.
    """
    namespace = namespace or Namespace.MAIN
    if name.endswith('-test'):
        stack = name[: -len('-test')]
        namespace = Namespace.TEST
    else:
        stack = name
    return stack, namespace


def build_dataset_name(stack: str, namespace: Namespace) -> str:
    """
    Dataset name is suffixed with "-test" for a test dataset (matching the
    corresponding sample-metadata project).
    """
    is_test = namespace != Namespace.MAIN
    return stack + ('-test' if is_test else '')


class Dataset(Target):
    """
    Represents a CPG dataset in a particular namespace: main or test.

    Each `dataset` at the CPG corresponds to
    * one GCP project: https://github.com/populationgenomics/team-docs/tree/main/storage_policies
    * one Pulumi stack: https://github.com/populationgenomics/analysis-runner/tree/main/stack
    * two sample metadata projects: main and test (the latter has a `-test` ending).

    An object of this class is parametrised by a dataset name and a namespace,
    meaning that it matches exactly one GCP project, exactly one stack, and exactly
    one sample metadata project.

    An object has two ID-like fields: `stack` and `name`:
    * `stack` is the name of the dataset (matches names of a GCP project or
       a Pulumi stack), e.g. "seqr", "hgdp".
    * `name` is the name of the namespace-specific sample-metadata project,
       e.g. "seqr", "seqr-test", "hgdp", "hgdp-test".
    """

    def __init__(
        self,
        name: str,
        namespace: Namespace | None = None,
        cohort: Cohort | None = None,
    ):
        super().__init__()
        self._sample_by_id: dict[str, Sample] = {}
        self.stack, self.namespace = parse_stack(name, namespace)
        self.cohort = cohort

    @staticmethod
    def create(
        name: str,
        namespace: Namespace,
    ) -> 'Dataset':
        """
        Create a dataset.
        """
        # Normalising the dataset's name:
        name = build_dataset_name(*parse_stack(name, namespace))
        return Dataset(
            name=name,
            namespace=namespace,
        )

    @property
    def is_test(self) -> bool:
        """
        If it's a test dataset.
        """
        return self.namespace != Namespace.MAIN

    @property
    def name(self) -> str:
        """
        Name is suffixed with "-test" for a test dataset (matching the
        corresponding sample-metadata project).
        """
        return build_dataset_name(self.stack, self.namespace)

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
        return (
            '' 
            if not self.cohort or self.cohort.sequencing_type == SequencingType.GENOME 
            else self.cohort.sequencing_type.value
        )

    def prefix(self, **kwargs) -> Path:
        """
        The primary storage path.
        """
        return to_path(dataset_path(
            self._seq_type_subdir(),
            dataset=self.stack,
            **kwargs,
        ))

    def tmp_prefix(self, **kwargs) -> Path:
        """
        Storage path for temporary files.
        """
        return to_path(dataset_path(
            self._seq_type_subdir(),
            dataset=self.stack,
            category='tmp',
            **kwargs,
        ))

    def web_prefix(self, **kwargs) -> Path:
        """
        Path for files served by an HTTP server Matches corresponding URLs returns by
        self.web_url() URLs.
        """
        return to_path(dataset_path(
            self._seq_type_subdir(),
            dataset=self.stack,
            category='web',
            **kwargs,
        ))

    def web_url(self, **kwargs) -> str | None:
        """
        URLs matching self.storage_web_path() files serverd by an HTTP server. 
        """
        return web_url(
            self._seq_type_subdir(),
            dataset=self.stack,
            **kwargs,
        )

    def add_sample(
        self,
        id: str,  # pylint: disable=redefined-builtin
        external_id: str | None = None,
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Optional['Sex'] = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input_by_seq_type: dict[SequencingType, AlignmentInput] | None = None,
    ) -> 'Sample':
        """
        Create a new sample and add it to the dataset.
        """
        if id in self._sample_by_id:
            logger.debug(f'Sample {id} already exists in the dataset {self.name}')
            return self._sample_by_id[id]

        s = Sample(
            id=id,
            dataset=self,
            external_id=external_id,
            participant_id=participant_id,
            meta=meta,
            sex=sex,
            pedigree=pedigree,
            alignment_input_by_seq_type=alignment_input_by_seq_type,
        )
        self._sample_by_id[id] = s
        return s

    def get_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Get dataset's samples. Include only "active" samples, unless only_active=False
        """
        return [
            s for sid, s in self._sample_by_id.items() if 
            (s.active or not only_active)
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

    def make_ped_file(self, tmp_bucket: Path | None = None) -> Path:
        """
        Create a PED file for all samples
        """
        datas = []
        for sample in self.get_samples():
            if sample.pedigree:
                datas.append(sample.pedigree.get_ped_dict())
        df = pd.DataFrame(datas)

        ped_path = (tmp_bucket or self.tmp_prefix()) / f'{self.name}.ped'
        with ped_path.open('w') as fp:
            df.to_csv(fp, sep='\t', index=False)

        return ped_path


class Sex(Enum):
    """
    Sex as in PED format
    """

    UNKNOWN = 0
    MALE = 1
    FEMALE = 2

    @staticmethod
    def parse(sex: str | None) -> 'Sex':
        """
        Parse a string into a Sex object.
        """
        if sex:
            if sex.lower() in ('m', 'male', '1'):
                return Sex.MALE
            if sex.lower() in ('f', 'female', '2'):
                return Sex.FEMALE
            if sex.lower() in ('u', 'unknown', '0'):
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
        alignment_input_by_seq_type: dict[SequencingType, AlignmentInput] | None = None,
    ):
        super().__init__()
        self.id = id
        self._external_id = external_id
        self.dataset = dataset
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.pedigree: PedigreeInfo | None = pedigree
        if sex:
            self.pedigree = PedigreeInfo(
                sample=self,
                fam_id=self.participant_id,
                sex=sex,
            )
        self.alignment_input_by_seq_type: dict[SequencingType, AlignmentInput] = \
            alignment_input_by_seq_type or dict()

    def __repr__(self):
        values = {
            'participant': self._participant_id if self._participant_id else '',
            'forced': str(self.forced),
            'active': str(self.active),
            'meta': str(self.meta),
            'alignment_inputs': ','.join([
                f'{seq_t.value}: {al_inp}' 
                for seq_t, al_inp in self.alignment_input_by_seq_type.items()
            ]),
            'pedigree': self.pedigree if self.pedigree else '',
        }
        retval = f'Sample({self.dataset.name}/{self.id}'
        if self._external_id:
            retval += f'|{self._external_id}'
        return retval + ''.join(f', {k}={v}' for k, v in values.items())

    def __str__(self):
        ai_tag = ''
        for seq_type, alignment_input in self.alignment_input_by_seq_type.items():
            ai_tag += f'|SEQ={seq_type.value}:'
            if isinstance(alignment_input, CramPath):
                if alignment_input.is_bam:
                    ai_tag += 'CRAM'
                else:
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

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sample, corresponding
        a PED file entry.
        """
        if self.pedigree:
            return self.pedigree.get_ped_dict(use_participant_id)
        return {
            'Family.ID': self.participant_id if use_participant_id else self.id,
            'Individual.ID': self.participant_id if use_participant_id else self.id,
            'Father.ID': '0',
            'Mother.ID': '0',
            'Sex': '0',
            'Phenotype': '0',
        }

    def get_cram_path(self) -> CramPath:
        """
        Path to a CRAM file. Not checking its existence here.
        """
        return CramPath(self.dataset.prefix() / 'cram' / f'{self.id}.cram')

    def get_gvcf_path(self) -> GvcfPath:
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
        return {
            'dataset': self.dataset.name, 
            'sample': self.id,
        }

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
    sex: Sex
    fam_id: str | None = None
    phenotype: str | None = None
    dad: Sample | None = None
    mom: Sample | None = None

    def get_ped_dict(self, use_participant_id: bool = False) -> dict:
        """
        Returns a dictionary of pedigree fields for this sample, corresponding
        a PED file entry.
        """

        def _get_id(_s: Sample | None):
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
            'Phenotype': self.phenotype,
        }
