#!/usr/bin/env python

# Libraries
from setuptools import setup, find_packages

# Description
with open('README.md') as fd:
    long_description = fd.read()
#end with

# Requirements
with open('requirements.txt') as fr:
    set_parsed = fr.read().splitlines()
#end with

# Set requires
install_requires = [req.strip() for req in set_parsed]
tests_requires = [
    'pytest'
]

# Setup
setup(
    name = 'granite-suite',
    version = open('granite/_version.py').readlines()[-1].split()[-1].strip("\"'"),
    description = 'granite is a collection of software to call, filter and work with genomic variants',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/dbmi-bgm/granite',
    classifiers=[
            'License :: OSI Approved :: MIT License',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3'
            ],
    install_requires=install_requires,
    setup_requires=install_requires,
    tests_require=tests_requires,
    include_package_data=True,
    packages = find_packages(),
    author = 'Michele Berselli',
    author_email = 'berselli.michele@gmail.com',
    license = 'MIT',
    python_requires = '>=3.6',
    entry_points = {
        'console_scripts': [
            'granite = granite.granite:main',
        ]
    }
    )
