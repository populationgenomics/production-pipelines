"""
Pipeline stage targets: cohort, dataset, sample.
"""

import logging
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from cpg_pipes.storage import Path

from cpg_pipes.pipeline.analysis import AlignmentInput, CramPath, GvcfPath
from cpg_pipes.storage import Namespace, StorageProvider
from cpg_pipes.pipeline.sequence import SmSequence
from cpg_pipes.pipeline.target import Target

logger = logging.getLogger(__file__)


class Cohort(Target):
    """
    Represents a "cohort" target - all samples from all datasets in the pipeline
    """
    def __init__(
        self, 
        name: str, 
        analysis_dataset_name: str,
        namespace: Namespace,
        storage_provider: StorageProvider | None = None,
    ):
        """
        @param name: name of the cohort
        @param analysis_dataset_name: deferring creation of the analysis
        dataset in case if it is also a part of the cohort datasets.
        """
        super().__init__()
        self.name = name
        self.storage_provider = storage_provider
        self.analysis_dataset = Dataset(
            name=analysis_dataset_name, 
            cohort=self,
            namespace=namespace,
            storage_provider=storage_provider or self.storage_provider,
        )
        self._datasets_by_name: dict[str, Dataset] = {}

    def __repr__(self):
        return self.name

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return f'Cohort("{self.name}", {len(self.get_datasets())} datasets)'

    def get_datasets(self, only_active: bool = True) -> list['Dataset']:
        """
        Gets list of all datasets. 
        Include only "active" datasets (unless only_active is False)
        """
        return [
            ds for k, ds in self._datasets_by_name.items() 
            if (ds.active or not only_active)
        ]
    
    def get_dataset_by_name(
        self, name: str, only_active: bool = True
    ) -> Optional['Dataset']:
        """
        Get dataset by name.
        Include only "active" datasets (unless only_active is False)
        """
        ds_by_name = {d.name: d for d in self.get_datasets(only_active)}
        return ds_by_name.get(name)

    def get_all_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Gets a flat list of all samples from all datasets.
        Include only "active" samples (unless only_active is False)
        """
        all_samples = []
        for proj in self.get_datasets(only_active=only_active):
            all_samples.extend(proj.get_samples(only_active))
        return all_samples

    def get_all_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Gets a flat list of CPG IDs for all samples from all datasets.
        """
        return [s.id for s in self.get_all_samples(only_active=only_active)]

    def add_dataset(
        self,
        name: str, 
        namespace: Namespace | None = None,
        storage_provider: StorageProvider | None = None,
    ) -> 'Dataset':
        """
        Create a dataset and add it into the cohort.
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
                cohort=self,
                namespace=namespace,
                storage_provider=storage_provider or self.storage_provider,
            )

        self._datasets_by_name[ds.name] = ds
        return ds


def parse_stack(
    name: str, 
    namespace: Namespace | None = None
) -> tuple[str, Namespace]:
    """
    Input `name` can be either e.g. "seqr" or "seqr-test". The latter will be 
    resolved to stack="seqr" and is_test=True, unless `namespace` is provided 
    explicitly.

    Returns the stack id and a corrected namespace.
    """
    namespace = namespace or Namespace.MAIN
    if name.endswith('-test'):
        stack = name[:-len('-test')]
        if namespace == Namespace.MAIN:
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
        cohort: Cohort,
        namespace: Namespace | None = None,
        storage_provider: StorageProvider | None = None,
    ):
        super().__init__()
        self.cohort = cohort
        
        self._storage_provider = storage_provider

        self._sample_by_id: dict[str, Sample] = {}
        
        self.stack, self.namespace = parse_stack(name, namespace)

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
        return f'Dataset("{self.name}")'
    
    def __repr__(self):
        return f'Dataset("{self.name}", {len(self.get_samples())} samples)'

    def __str__(self):
        return f'{self.name} ({len(self.get_samples())} samples)'

    @property
    def storage_provider(self) -> StorageProvider:
        """
        Storage provider required to get bucket paths for the dataset.
        """
        if not self._storage_provider:
            raise ValueError(
                '_storage_provider is not set. storage_provider must be passed to the '
                'Dataset() constructor before calling Dataset.get_bucket()'
            )
        return self._storage_provider

    def get_bucket(self, **kwargs) -> Path:
        """
        The primary dataset bucket (-main or -test).
        """
        return self.storage_provider.get_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs,
        )

    def get_tmp_bucket(self, **kwargs) -> Path:
        """
        The tmp bucket (-main-tmp or -test-tmp)
        """
        return self.storage_provider.get_tmp_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def get_analysis_bucket(self, **kwargs) -> Path:
        """
        Get analysis bucket (-main-analysis or -test-analysis)
        """
        return self.storage_provider.get_analysis_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def get_web_bucket(self, **kwargs) -> Path:
        """
        Get web bucket (-main-web or -test-web)
        """
        return self.storage_provider.get_web_bucket(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )

    def get_web_url(self, **kwargs) -> str | None:
        """
        Get web base URL.
        """
        return self.storage_provider.get_web_url(
            dataset=self.stack, namespace=self.namespace, **kwargs
        )
    
    def add_sample(
        self, 
        id: str,  # pylint: disable=redefined-builtin
        external_id: str, 
        participant_id: str | None = None,
        seq: SmSequence | None = None,
        sex: Optional['Sex'] = None,
        pedigree: Optional['PedigreeInfo'] = None,
        **kwargs
    ) -> 'Sample':
        """
        Create a new sample and add it into the dataset.
        """
        if id in self._sample_by_id:
            logger.debug(f'Sample {id} already exists in the dataset {self.name}')
            return self._sample_by_id[id]

        s = Sample(
            id=id, 
            external_id=external_id,
            participant_id=participant_id,
            seq=seq,
            sex=sex,
            pedigree=pedigree,
            dataset=self,
            meta=kwargs,
        )
        self._sample_by_id[id] = s
        return s
    
    def get_samples(self, only_active: bool = True) -> list['Sample']:
        """
        Get dataset's samples. Inlcude only "active" samples, unless only_active=False
        """
        return [
            s for sid, s in self._sample_by_id.items() 
            if (s.active or not only_active)
        ]

    def get_sample_ids(self, only_active: bool = True) -> list[str]:
        """
        Get dataset's sample IDs. Inlcude only "active" samples, 
        unless only_active=False.
        """
        return [s.id for s in self.get_samples(only_active=only_active)]


class Sample(Target):
    """
    Represents a Sample.
    """

    def __init__(
        self,
        id: str,  # pylint: disable=redefined-builtin
        external_id: str,
        dataset: 'Dataset',  # type: ignore  # noqa: F821
        participant_id: str | None = None,
        meta: dict | None = None,
        sex: Optional['Sex'] = None,
        seq: SmSequence | None = None,
        pedigree: Optional['PedigreeInfo'] = None,
        alignment_input: AlignmentInput | None = None
    ):
        super().__init__()
        self.id = id
        self.external_id = external_id
        self.dataset = dataset
        self._participant_id = participant_id
        self.meta: dict = meta or dict()
        self.seq = seq
        self.pedigree: PedigreeInfo | None = pedigree
        if sex:
            self.pedigree = PedigreeInfo(
                sample=self,
                fam_id=self.participant_id,
                sex=sex,
            )
        self._alignment_input = alignment_input

    def __repr__(self):
        return (
            f'Sample({self.dataset.name}/{self.id}|{self.external_id}' +
            (f', participant={self._participant_id}'
             if self._participant_id else '') +
            f', forced={self.forced}' +
            f', active={self.active}' +
            f', meta={self.meta}' +
            (f', seq={self.seq}' if self.seq else '') +
            (f', alignment_input={self._alignment_input}'
             if self._alignment_input else '') +
            (f', pedigree={self.pedigree}' if self.pedigree else '') +
            f')'
        )

    def __str__(self):
        ai_tag = ''
        if self._alignment_input:
            if isinstance(self._alignment_input, CramPath):
                if self._alignment_input.is_bam:
                    ai_tag = f'|SEQ=CRAM'
                else:
                    ai_tag = f'|SEQ=BAM'
            else:
                ai_tag = f'|SEQ={len(self._alignment_input)}FQS'

        return f'Sample({self.dataset.name}/{self.id}|{self.external_id}{ai_tag})'

    @property
    def alignment_input(self) -> AlignmentInput | None:
        """
        Returns input for (re-)alignment.
        """
        if self._alignment_input:
            return self._alignment_input

        if self.seq:
            return self.seq.alignment_input
        return None

    @alignment_input.setter
    def alignment_input(self, value: AlignmentInput):
        """
        Set alignment_input
        """
        self._alignment_input = value

    @property
    def participant_id(self) -> str:
        """
        Get participant's ID. 
        Uses external_id whenever participant_id is not available in the DB
        """
        return self._participant_id or self.external_id

    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return f'Sample("{self.id})"'

    def get_ped_dict(self, use_participant_id: bool = False) -> dict[str, str]:
        """
        Returns a dictionary of pedigree fields for this sample, corresponging
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
        return CramPath(self.dataset.get_bucket() / 'cram' / f'{self.id}.cram')

    def get_gvcf_path(self) -> GvcfPath:
        """
        Path to a GVCF file. Not checking its existence here.
        """
        return GvcfPath(self.dataset.get_bucket() / 'gvcf' / f'{self.id}.g.vcf.gz')


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
        return Sex.UNKNOWN


@dataclass
class PedigreeInfo:
    """
    Pedigree relationsips with other samples in the cohort, and other PED data
    """
    sample: Sample
    sex: Sex
    fam_id: str | None = None
    phenotype: str | None = None
    dad: Sample | None = None
    mom: Sample | None = None

    def get_ped_dict(self, use_participant_id: bool = False) -> dict:
        """
        Returns a dictionary of pedigree fields for this sample, corresponging
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


class Pair(Target):
    """
    Pair of samples
    """
    @property
    def target_id(self) -> str:
        """Unique target ID"""
        return f'{self.s1.target_id}:{self.s2.target_id}'

    def __init__(self, s1: Sample, s2: Sample):
        super().__init__()
        self.s1 = s1
        self.s2 = s2
