#!/usr/bin/env python


from setuptools import setup, find_packages


requires = [
]


setup(
    name="granite",
    version=open("granite/_version.py").readlines()[-1].split()[-1].strip("\"'"),
    description='granite is a collection of software to call, filter and work with genomic variants',
    entry_points = {
        'console_scripts': [
            'granite = granite.granite:main',
        ]
    },
    packages=find_packages(),
    author='Michele Berselli',
    author_email='berselli.michele@gmail.com'
    )
