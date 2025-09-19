#!/usr/bin/env python3

import os
import sys
import time
import shutil
import argparse
import numpy as np

from ase.lattice import FCC
from ase.dft.kpoints import *

parser = argparse.ArgumentParser()
parser.add_argument('--clean', action='store_true')
parser.add_argument('--print_path', action='store_true')
args = parser.parse_args()

def clean(dir_output=os.getcwd(), keep=['Pt.poscar', 'Pt.sh', 'Pt_band_.abi', 'Pt_crpa_.abi', 'log', 'run.py']):
	for fn in set(os.listdir(dir_output)) ^ set(keep):
		fn = f'{dir_output}/{fn}'
		if os.path.isfile(fn): os.remove(fn)
		else: shutil.rmtree(fn)

def print_path(path='GXWKGLUWLKX'):
	coords = get_special_points(FCC(3).tocell())
	for kpt in path:
		print(*coords[kpt], f'#{kpt}')

if args.clean: clean()
elif args.print_path: print_path()
else:
	sys.exit(1)
	parser.print_help()
