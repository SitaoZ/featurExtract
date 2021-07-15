# -*- coding: utf-8 -*-
import os

def readme():
    with open("README.rst") as f:
        return f.read()

from setuptools import setup 
from featurExtract.version import __version__

setup(
    name='featurExtract',
    version=__version__,
    keywords='genome feature, extract',
    description='Extract genome ferature sequence for biologists',
    long_description_content_type='text/markdown',
    long_description=readme(),
    entry_points = {'console_scripts': [
                       'featurExtract=featurExtract.command_gff:main',
                       'genBankExtract=featurExtract.command_gb:main'
                   ]},
    author='zhusitao',
    author_email='zhusitao1990@163.com',
    url='https://github.com/SitaoZ/featurExtract.git',
    include_package_data=True,
    packages=['featurExtract'],
    license='MIT',
    install_requires = ['argparse>=1.1', 
                        'pandas>=1.2.4', 
                        'gffutils>=0.10.1',
                        'setuptools>=49.2.0',
                        'biopython>=1.78'])
