"""
Stage that generates a CRAM file.
"""

import logging

from cpg_utils.config import get_config

from .. import Path
from ..jobs.align import Aligner, MarkDupTool
from ..jobs.cram_qc import samtools_stats, picard_wgs_metrics, verify_bamid
from ..targets import Sample
from ..pipeline import stage, SampleStage, StageInput, StageOutput
from ..jobs import align, somalier
from ..types import CramPath

logger = logging.getLogger(__file__)


@stage(analysis_type='cram')
class Align(SampleStage):
    """
    Align or re-align input data to produce a CRAM file
    """

    def expected_outputs(self, sample: Sample) -> dict[str, Path]:
        """
        Stage is expected to generate a CRAM file and a corresponding index.
        """
        qc_dict = {
            key: sample.dataset.prefix() / 'qc' / key / f'{sample.id}{suf}'
            for key, suf in {
                'markduplicates_metrics': '.markduplicates-metrics',
                'samtools_stats': '_samtools_stats.txt',
                'picard_wgs_metrics': '_picard_wgs_metrics.csv',
                'verify_bamid': '_verify_bamid.selfSM',
            }.items()
        }
        return {
            'cram': sample.get_cram_path().path,
            'somalier': sample.get_cram_path().somalier_path,
        } | qc_dict

    @staticmethod
    def _get_cram_reference_from_version(cram_version) -> str:
        """
        Get the reference used for the specific cram_version,
        so that bazam is able to correctly decompress the reads
        """
        cram_version_map = get_config()['workflow'].get('cram_version_reference', {})
        if cram_version in cram_version_map:
            return cram_version_map[cram_version]
        raise ValueError(
            f'Unrecognised cram_version: "{cram_version}", expected one of: {", ".join(cram_version_map.keys())}'
        )

    def queue_jobs(self, sample: Sample, inputs: StageInput) -> StageOutput | None:
        """
        Using the "align" function implemented in the `jobs` module.
        Checks the `realign_from_cram_version` pipeline config argument, and
        prioritises realignment from CRAM vs alignment from FASTQ if it's set.
        """
        alignment_input = sample.alignment_input_by_seq_type.get(
            self.cohort.sequencing_type
        )
        if realign_cram_ver := get_config()['workflow'].get(
            'realign_from_cram_version'
        ):
            if (
                path := (
                    sample.dataset.prefix()
                    / 'cram'
                    / realign_cram_ver
                    / f'{sample.id}.cram'
                )
            ).exists():
                logger.info(f'Realigning from {realign_cram_ver} CRAM {path}')
                alignment_input = CramPath(
                    path,
                    reference_assembly=self._get_cram_reference_from_version(
                        realign_cram_ver
                    ),
                )

        if alignment_input is None or (
            get_config()['workflow'].get('check_inputs')
            and not alignment_input.exists()
        ):
            if get_config()['workflow'].get('skip_samples_with_missing_input'):
                logger.error(f'No alignment inputs, skipping sample {sample}')
                sample.active = False
                return self.make_outputs(sample, skipped=True)  # return empty output
            else:
                return self.make_outputs(
                    target=sample, error_msg=f'No alignment input found'
                )
        assert alignment_input

        jobs = align.align(
            b=self.b,
            alignment_input=alignment_input,
            output_path=sample.get_cram_path(),
            out_markdup_metrics_path=self.expected_outputs(sample)[
                'markduplicates_metrics'
            ],
            sample_name=sample.id,
            job_attrs=self.get_job_attrs(sample),
            overwrite=not get_config()['workflow'].get('check_intermediates'),
            realignment_shards_num=get_config()['workflow'].get(
                'realignment_shards_num', align.DEFAULT_REALIGNMENT_SHARD_NUM
            ),
            aligner=Aligner.DRAGMAP,
            markdup_tool=MarkDupTool.PICARD,
            sequencing_type=self.cohort.sequencing_type,
        )

        for key, qc_func in {
            'samtools_stats': samtools_stats,
            'picard_wgs_metrics': picard_wgs_metrics,
            'verify_bamid': verify_bamid,
            'somalier': somalier.extact_job,
        }:
            j = qc_func(
                self.b,
                sample.get_cram_path(),
                job_attrs=self.get_job_attrs(sample),
                output_path=self.expected_outputs(sample)[key],
                sequencing_type=self.cohort.sequencing_type,
                overwrite=not get_config()['workflow'].get('check_intermediates'),
            )
            if j:
                jobs.append(j)

        return self.make_outputs(sample, data=self.expected_outputs(sample), jobs=jobs)
