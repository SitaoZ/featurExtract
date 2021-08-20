# -*- coding: utf-8 -*-
import os

def readme():
    with open("README.md") as f:
        long_description = f.read()
        return long_description

from setuptools import setup 
from featurExtract.version import __version__

PACKAGES = [
    "featurExtract",
    "featurExtract.commands",
    "featurExtract.database",
    "featurExtract.utils"
]
setup(
    name='featurExtract',
    version=__version__,
    keywords='genome feature, extract',
    description='Extract genome ferature sequence for biologists',
    long_description=readme(),
    long_description_content_type='text/markdown',
    entry_points = {'console_scripts': [
                       'featurExtract=featurExtract.command_gff:main',
                       'genBankExtract=featurExtract.command_gb:main'
                   ]},
    author='zhusitao',
    author_email='zhusitao1990@163.com',
    url='https://github.com/SitaoZ/featurExtract.git',
    include_package_data=True, # done via MANIFEST.in under setuptools
    packages=PACKAGES,
    license='MIT',
    install_requires = ['argparse>=1.1', 
                        'pandas>=1.2.4', 
                        'gffutils>=0.10.1',
                        'setuptools>=49.2.0',
                        'biopython>=1.78'],
    python_requires=">=3.7.6")
