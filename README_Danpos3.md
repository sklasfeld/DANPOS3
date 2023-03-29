# DANPOS3

For questions on this code, please feel free to checkout the DANPOS project forum at:
https://groups.google.com/forum/#!forum/danpos

## Background

This repo contains updated code from DANPOS2.
* This package was originally forked from code from danpos-2.2.2 (Jun 15 10:20) which is dependent on python 2.7
 * https://sites.google.com/site/danposdoc/download/danpos-2.2.2.tgz
* For information on how to use Danpos3 commands please see the Danpos2 manual at :
https://sites.google.com/site/danposdoc/ where this code originated from.

## Summary

This code was written to update umap to work with Python3 and fix other bugs.

## Installation
To download the repository
```
git clone https://github.com/sklasfeld/DANPOS3.git
```

To install the exact versions of the dependencies use type:
```
cd DANPOS3
pip install -r requirements.txt
```

[To be able to execute these scripts from other directories
easily, please set your current working directory to your
$PATH variable.](https://opensource.com/article/17/6/set-path-linux)

### Package and Library versions
* Python 3.7.6
* R version 4.0.1
* samtools 1.7 using htslib 1.7
* Python Libraries
  * rpy2 3.3.3
  * argparse 1.1
  * numpy 1.18.5
  * pysam 0.16.0.1
  * scipy 1.5.4

To test your environment for these packages, use the script `print_versions.py`
