#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='cpg-pipes',
    version='0.2.11',
    description=(
        'Hail Batch pipelines for large cohort and rare deseases projects'  
    ),
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=['cpg_pipes'],
    package_data={'cpg_pipes': ['filter_cutoffs.yaml']},
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', fname) 
        for fname in os.listdir('scripts') if fname.endswith('.py')
    ],
    install_requires=[
        'click',
        'pandas',
        'cpg-hail',
        'cpg-gnomad',   # github.com/populationgenomics/gnomad_methods
        'google-cloud-storage',
        'fsspec',
        'sample-metadata',
        'analysis-runner',
    ],
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
