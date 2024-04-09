#!/usr/bin/env python

from setuptools import find_packages, setup


def read_reqs(filename: str) -> list[str]:
    """
    Read requirements from a file, return as a list
    Args:
        filename (str): the requirements file to parse

    Returns:
        list[str]: the requirements
    """
    with open(filename, encoding='utf-8') as filehandler:
        return [line.strip() for line in filehandler if line.strip() and not line.startswith('#')]


setup(
    name='cpg-workflows',
    # This tag is automatically updated by bumpversion
    version='1.22.7',
    description='CPG workflows for Hail Batch',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/production-pipelines',
    license='MIT',
    packages=find_packages(),
    install_requires=read_reqs('requirements.txt'),
    extras_require={'dev': read_reqs('requirements-dev.txt')},
    package_data={'cpg_workflows': ['defaults.toml']},
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
