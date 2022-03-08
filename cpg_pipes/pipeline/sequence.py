"""
Sample-metadata DB Sequence entry.
"""

import logging

from cpg_pipes import buckets
from cpg_pipes.pipeline.analysis import AlignmentInput

logger = logging.getLogger(__file__)
logging.basicConfig(format='%(levelname)s (%(name)s %(lineno)s): %(message)s')
logger.setLevel(logging.INFO)


class SmSequence:
    """
    Sample-metadata DB "Sequence" entry.

    See sample-metadata for more details: https://github.com/populationgenomics/sample-metadata
    """
    def __init__(
        self, 
        id: str,
        sample_id: str, 
        meta: dict,
        alignment_input: AlignmentInput|None
    ):
        self.id = id
        self.sample_id = sample_id
        self.meta = meta
        self.alignment_input = alignment_input

    @staticmethod
    def parse(data: dict, check_existence: bool) -> 'SmSequence':
        req_keys = ['id', 'sample_id', 'meta']
        if any(k not in data for k in req_keys):
            for key in req_keys:
                if key not in data:
                    logger.error(f'"Sequence" data does not have {key}: {data}')
            raise ValueError(f'Cannot parse SMDB Sequence {data}')

        return SmSequence(
            id=data['id'], 
            sample_id=data['sample_id'], 
            meta=data['meta'], 
            alignment_input=SmSequence._parse_reads(
                meta=data['meta'], 
                check_existence=check_existence
            ),
        )

    @staticmethod
    def _parse_reads(  # pylint: disable=too-many-return-statements
        meta: dict,
        check_existence: bool,
    ) -> AlignmentInput|None:
        """
        Parse a AlignmentInput object from the meta dictionary.

        :param check_existence: check if fastq/crams exist on buckets. 
        Default value is pulled from self.smdb and can be overwridden.
        """
        reads_data = meta.get('reads')
        reads_type = meta.get('reads_type')
        if not reads_data:
            logger.error(f'No "meta/reads" field in meta')
            return None
        if not reads_type:
            logger.error(f'No "meta/reads_type" field in meta')
            return None
        supported_types = ('fastq', 'bam', 'cram')
        if reads_type not in supported_types:
            logger.error(f'ERROR: "reads_type" is expected to be one of {supported_types}')
            return None
    
        if reads_type in ('bam', 'cram'):
            if len(reads_data) > 1:
                logger.error('Supporting only single bam/cram input')
                return None
    
            bam_path = reads_data[0]['location']
            if not (bam_path.endswith('.cram') or bam_path.endswith('.bam')):
                logger.error(
                    f'ERROR: expected the file to have an extention .cram or .bam,'
                    f'got: {bam_path}'
                )
                return None
            if check_existence and not buckets.file_exists(bam_path):
                logger.error(f'ERROR: index file doesn\'t exist: {bam_path}')
                return None
    
            # Index:
            if not reads_data[0].get('secondaryFiles'):
                logger.error(
                    f'ERROR: bam/cram input is expected to have '
                    f'a non-empty list field "secondaryFile" section with indices'
                )
                return None
            index_path = reads_data[0]['secondaryFiles'][0]['location']
            if (
                bam_path.endswith('.cram')
                and not index_path.endswith('.crai')
                or bam_path.endswith('.bai')
                and not index_path.endswith('.bai')
            ):
                logger.error(
                    f'ERROR: expected the index file to have an extention '
                    f'.crai or .bai, got: {index_path}'
                )
            if check_existence and not buckets.file_exists(index_path):
                logger.error(f'ERROR: index file doesn\'t exist: {index_path}')
                return None
    
            return AlignmentInput(cram_path=bam_path, index_path=index_path)

        else:
            fqs1 = []
            fqs2 = []
            for lane_data in reads_data:
                assert len(lane_data) == 2, lane_data
                if check_existence and not buckets.file_exists(lane_data[0]['location']):
                    logger.error(
                        f'ERROR: read 1 file doesn\'t exist: {lane_data[0]["location"]}'
                    )
                    return None
                if check_existence and not buckets.file_exists(lane_data[1]['location']):
                    logger.error(
                        f'ERROR: read 2 file doesn\'t exist: {lane_data[1]["location"]}'
                    )
                    return None
    
                fqs1.append(lane_data[0]['location'])
                fqs2.append(lane_data[1]['location'])
            return AlignmentInput(fqs1=fqs1, fqs2=fqs2)
