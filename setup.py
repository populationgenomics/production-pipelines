#!/usr/bin/env python

import os
from os.path import join
import setuptools


setuptools.setup(
    name='cpg-pipes',
    version='0.3.5',
    description='Hail Batch bioinformatics pipelines',
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
        for fname in os.listdir('scripts')
        if fname.endswith('.py')
    ],
    install_requires=[
        'cpg-utils',
        'pandas',
        'hail>=0.2.91',
        'cpg-gnomad',  # github.com/populationgenomics/gnomad_methods
        'peddy',
        'google-cloud-storage',
        'google-cloud-secret-manager',
        'fsspec',
        'sample-metadata',
        'analysis-runner',
        'cloudpathlib[all]',
        'coloredlogs',
        'types-PyYAML',  # https://mypy.readthedocs.io/en/stable/getting_started.html#library-stubs-and-typeshed
        'python-slugify',
    ],
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
