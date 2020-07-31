import sys
import pkg_resources
import os

print("\nPython version")
print (sys.version)
required = {'rpy2','argparse','numpy','pysam'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed
available = [pkg for pkg in required if pkg not in missing]
print("-------------------------------------------------\n")

print(os.system("R --version"))
print("-------------------------------------------------\n")

print(os.system("samtools --version"))
print("-------------------------------------------------\n")

print("AVAILABLE python packages:")
for pkg in available:
	if pkg == "rpy2":
		import rpy2
		print(pkg + ": "+ str(rpy2.__version__))
	elif pkg == "argparse":
		import argparse
		print (pkg + ": "+ str(argparse.__version__))
	elif pkg == "numpy":
		import numpy as np
		print("numpy: "+ str(np.__version__))
	else:
		import pysam
		print("pysam: "+ str(pysam.__version__))

if len(missing) > 0:
	print("\nMISSING python packages:")
	print(missing)
else:
	print("\nGreat job! No missing python packages!")
