#!/usr/bin/env python3

"""
Entry point to run the pipeline.
"""

import coloredlogs
from cpg_utils.workflows.workflow import get_workflow
from stages.seqr_loader import MtToEs

fmt = '%(asctime)s %(levelname)s (%(name)s %(lineno)s): %(message)s'
coloredlogs.install(level='DEBUG', fmt=fmt)


def main():
    get_workflow().run(stages=[MtToEs])


if __name__ == '__main__':
    main()
