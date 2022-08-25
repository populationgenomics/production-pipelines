#!/usr/bin/env python

import setuptools


setuptools.setup(
    name='cpg-seqr-loader',
    version='0.3.10',
    description='CPG Seqr loading pipeline',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=['jobs', 'stages']
    + ['jobs.' + p for p in sorted(setuptools.find_packages('./jobs'))]
    + ['stages.' + p for p in sorted(setuptools.find_packages('./stages'))],
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['main=main'],
        'stages': [
            'MtToEs = stages.seqr_loader:MtToEs',
            'CramMultiQC = stages.multiqc:CramMultiQC',
            'GvcfMultiQC = stages.multiqc:GvcfMultiQC',
        ],
    },
    install_requires=[
        'cpg-utils',
        'sample-metadata',
        'analysis-runner',
        'cpg-gnomad',  # https://github.com/populationgenomics/gnomad_methods
        'pytest',
        'pandas',
        'hail',
        'peddy',
        'google-cloud-storage',
        'google-cloud-secret-manager',
        'fsspec',
        'cloudpathlib[all]',
        'coloredlogs',
        'slack_sdk',
        'elasticsearch==8.*',
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
