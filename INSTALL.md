# Install Guide for DANPOS3

Please check the following instructions to complete your installation

## Prerequisites
The following are the recommended versions, but any version above
the following should also work.
* Python 3.7.6
* R version 4.0.1
* samtools 1.7 using htslib 1.7
* Python Libraries
  * rpy2 3.3.3
  * argparse 1.1
  * numpy 1.18.5
  * pysam 0.16.0.1

## Current environment

To see which versions are currently in your environment, use the command:
```
python scripts/print_versions.py
```

## Install from source

1. Download the repository
```
git clone https://github.com/sklasfeld/DANPOS3.git
```

2. Install the exact versions of the dependencies type:

```
cd DANPOS3
pip install -r requirements.txt
```

3. run the install script

To install the script globally:
```
python setup.py install
```
If you need to provide a nonstandard install prefix, or any other nonstandard options, you can provide many command line options to the install script. Use the –help option to see a brief list of available options:
```
python setup.py --help
```
For example, to install DANPOS3 locally under the HOME directory (/home/sklasfeld), use this command:
```
python setup.py install --prefix /home/sklasfeld
```

## Environmental Variables
If you locally installed DANPOS3, you may have to set up your `PYTHONPATH` and `PATH` environmental variables. The process for setting your environmental variables varies on each platform, but the general concept is the same. [See how to set your PATH on Linux.](https://opensource.com/article/17/6/set-path-linux)

### PYTHONPATH

To set up your PYTHONPATH environment variable, you'll need to add the value PREFIX/lib/pythonX.Y/site-packages to your existing PYTHONPATH. In this value, X.Y stands for the major–minor version of Python you are using (such as 3.7; you can find this using the `print_versions.py` script in the `scripts/` directory). PREFIX is the install prefix where you installed MACS. For example, if we installed the script to the HOME directory (/home/sklasfeld) using python 3.7, then we could export PYTHONPATH as such:
```
export PYTHONPATH=/home/sklasfeld/lib/python3.7/site-packages:$PYTHONPATH
```

Using Windows, you need to open up the system properties dialog and locate the tab labeled Environment. Add your value to the PYTHONPATH variable, or create a new PYTHONPATH variable if there isn't one already.

### PATH

Similar to how we exported the PYTHONPATH, to set up your PATH environment variable, you'll need to add the value PREFIX/bin to your existing PATH. For example:
```
export PATH=/home/sklasfeld/bin:$PATH
```