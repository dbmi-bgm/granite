#!/usr/bin/env python


from setuptools import setup, find_packages


requires = [
]


setup(
    name="granite",
    description='granite is a collection of software to call, filter and work with genomic variants',
    entry_points = {
        'console_scripts': [
            'granite = granite.granite:main',
        ]
    },
    packages=find_packages()
    )
