#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='cpg-pipes',
    version='0.2.0',
    description=(
        'Hail Batch pipelines for population genomics and rare deseases projects'  
    ),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=['cpg_production_pipelines'],
    package_data={'cpg_production_pipelines': ['filter_cutoffs.yaml']},
    include_package_data=True,
    zip_safe=False,
    scripts=[join('batches', fp) for fp in os.listdir('batches') if fp.endswith('.py')],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
