#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.30.3',
    description='CPG workflows for Hail Batch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'cpg-utils>=5.1.1',
        'cyvcf2==0.30.18',
        'analysis-runner>=2.43.3',
        'hail==0.2.133',  # Pin Hail at CPG's installed version
        'networkx>=2.8.3',
        'obonet>=0.3.1',  # for HPO parsing
        'grpcio-status>=1.48,<1.50',  # Avoid dependency resolution backtracking
        'onnx',
        'onnxruntime',
        'skl2onnx',
        'metamist>=6.9.0',
        'pandas',
        'peddy>=0.4.8',  # Avoid 0.4.7, which is incompatible
        'pyfaidx>=0.8.1.1',
        'fsspec',
        'slack_sdk',
        'elasticsearch==8.*',
        'coloredlogs',
        'bokeh',
        'numpy',
        'click',
        'tenacity',
        'toml',
    ],
    extras_require={
        'test': [
            'pytest',
            'pytest-xdist',
            'pytest-mock',
            'coverage',
        ],
    },
    package_data={
        'cpg_workflows': ['defaults.toml'],
    },
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
    entry_points={
        'console_scripts': [
            # script for modifying the content of a sniffles VCF, used in Long-Read SV pipeline
            'modify_sniffles = cpg_workflows.scripts.long_read_sniffles_vcf_modifier:cli_main',
            # for use in translating a MatrixTable to an ES index, first localising the MT
            'mt_to_es = cpg_workflows.scripts.mt_to_es_without_dataproc:main',
            # Generate new intervals from a MatrixTable
            'new_intervals_from_mt = cpg_workflows.scripts.generate_new_intervals:cli_main',
            'seqr_loader_cnv = cpg_workflows.scripts.seqr_loader_cnv:cli_main',
        ],
    },
)
