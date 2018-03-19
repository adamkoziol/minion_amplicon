#!/usr/bin/env python
from setuptools import setup
__author__ = 'adamkoziol'
setup(
    name="MinION Amplicon Pipeline",
    version="0.0.1",
    include_package_data=True,
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='Workflow for the automated binning and mapping of MinION reads to a set of targets',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md').read()
)
