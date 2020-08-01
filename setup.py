#!/usr/bin/env python3
# Time-stamp: 07-31-2020 14:45PM
# sklasfeld

"""Description: 
Setup script for DANPOS3: 
	A toolkit for Dynamic Analysis of Nucleosome 
	and Protein Occupancy by Sequencing, version 3
"""

from setuptools import setup

DESCRIPTION = ("DANPOS3: A toolkit for Dynamic Analysis of Nucleosome " +
  "and Protein Occupancy by Sequencing, version 3")

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='DANPOS3',
      version='3.1.1',
      description=DESCRIPTION,
      long_description=readme(),
      classifiers=[
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3"],
      url='https://github.com/sklasfeld/DANPOS3',
      author="Samantha Klasfeld",
      author_email='sjk314@gmail.com',
      license='GPLv3',
      packages=['danpos'],
      install_requires=[
          "argparse",
          "numpy",
          "pysam",
          "rpy2"],
      include_package_data=True,
      zip_safe=False)
