#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.25.5',
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
        'hail==0.2.130',  # Pin Hail at CPG's installed version
        'networkx>=2.8.3',
        'obonet>=0.3.1',  # for HPO parsing
        'grpcio-status>=1.62',  # Avoid dependency resolution backtracking
        'onnx',
        'onnxruntime',
        'skl2onnx',
        'metamist>=6.9.0',
        'pandas',
        'peddy',
        'fsspec',
        'slack_sdk',
        'elasticsearch==8.*',
        'coloredlogs',
        'bokeh',
        'numpy',
        'click',
        'toml',
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
