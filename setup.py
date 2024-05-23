#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.24.9',
    description='CPG workflows for Hail Batch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'cpg-utils>=5.0.4',
        'cyvcf2==0.30.18',
        'analysis-runner>=2.43.3',
        'hail==0.2.126',  # Temporarily pin a pre-job-groups Hail version
        'networkx>=2.8.3',
        'obonet>=0.3.1',  # for HPO parsing
        'grpcio-status>=1.48,<1.50',  # Avoid dependency resolution backtracking
        'onnx',
        'onnxruntime',
        'skl2onnx',
        'metamist>=6.9.0',
        'pandas',
        'peddy>=0.4.8',  # Avoid 0.4.7, which is incompatible
        'fsspec',
        'slack_sdk',
        'elasticsearch==8.*',
        'coloredlogs',
        'bokeh',
        'numpy',
        'click',
    ],
    extras_require={
        'test': [
            'pytest',
            'pytest-mock',
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
)
