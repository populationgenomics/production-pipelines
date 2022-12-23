#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.10.0',
    description='CPG workflows for Hail Batch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'cpg-utils',
        'analysis-runner>=2.36.4',
        'hail',
        'networkx',
        'sample-metadata>=5.0.1',
        'pandas',
        'peddy',
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
