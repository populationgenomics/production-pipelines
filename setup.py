#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.1.2',
    description='CPG workflows for Hail Batch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'cpg-utils',
        'networkx',
        'sample-metadata>=5.0.1',
        'gnomad',
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
        # Putting "hail" into an extra because by default we want to avoid
        # overriding the driver image CPG-fork based package:
        # https://github.com/populationgenomics/analysis-runner/blob/main/driver/Dockerfile.hail#L48-L57
        # "analysis-runner" also depends on "hail"
        'full': [
            'analysis-runner',
            'hail',
        ],
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
