#!/usr/bin/env python3

"""
Pipeline to align or realign available sequencing data to produce CRAMs
"""

from cpg_pipes.pipeline import pipeline_entry_point
from cpg_pipes.pipeline.pipeline import Pipeline
from cpg_pipes.stages.align import Align


@pipeline_entry_point(name='alignment')
def main(pipeline: Pipeline):
    """
    Pipeline to align or realign sequencing data to produce CRAMs
    """
    pipeline.run(stages=[Align])


if __name__ == '__main__':
    main()  # pylint: disable=E1120
