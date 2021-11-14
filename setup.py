#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='cpg-pipes',
    version='0.2.3',
    description=(
        'Hail Batch pipelines for population genomics and rare deseases projects'  
    ),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    # TODO: Can remove hail_scripts, lib, lib.model and package_dir after the oficial 
    # Hail package is updated to have https://github.com/hail-is/hail/pull/10863
    packages=[
        'cpg_pipes', 
        'hail_scripts', 
        'lib', 
        'lib.model'
    ],
    package_dir={
        'hail_scripts': 'hail-elasticsearch-pipelines/hail_scripts',
        'lib': 'hail-elasticsearch-pipelines/luigi_pipeline/lib',
        'lib.model': 'hail-elasticsearch-pipelines/luigi_pipeline/lib/model',
    },
    package_data={'cpg_pipes': ['filter_cutoffs.yaml']},
    include_package_data=True,
    zip_safe=False,
    scripts=[join('scripts', fname) for fname in os.listdir('scripts') if fname.endswith('.py')],
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
