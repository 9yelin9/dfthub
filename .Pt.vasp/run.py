#!/usr/bin/env python3

import os
import sys
import time
import shutil
import argparse
import numpy as np

from ase import io
from ase.calculators.vasp import Vasp
from ase.dft.kpoints import monkhorst_pack

parser = argparse.ArgumentParser()
parser.add_argument('--clean', action='store_true')
parser.add_argument('--wann', action='store_true')
parser.add_argument('--crpa', action='store_true')
args = parser.parse_args()

def clean(dir_output=os.getcwd(), keep=['Pt.poscar', 'Pt.sh', 'log', 'run.py']):
	for fn in set(os.listdir(dir_output)) ^ set(keep):
		fn = f'{dir_output}/{fn}'
		if os.path.isfile(fn): os.remove(fn)
		else: shutil.rmtree(fn)

def run_vasp(dir_output=os.getcwd(), np=4):
	atoms = io.read('Pt.poscar')

	atoms.calc = Vasp(
		directory=dir_output,
		command=f'mpirun -n {np} vasp_std',
		encut=300,  # Cutoff energy
		ismear=0,   # Gaussian smearing
		sigma=0.20, # Smearing width
		ediff=1e-6, # Convergence criteria for electronic steps
		icharg=2,   # Initial charge density
		nelm=100,   # Maximum iteration
		nsw=0,		# Ionic movement
		kpts=(4, 4, 4),
		lwave='.FALSE.',
		lcharg='.FALSE.',
	)

	t0 = time.time()
	atoms.get_potential_energy()
	tm, ts = divmod(int(time.time() - t0), 60)

	print(f'Done: {tm}m {ts}s', end='\n\n')

def run_wann(dir_output=os.getcwd()):
	pass

def run_crpa(dir_output=os.getcwd()):
	pass

if args.clean: clean()
elif args.wann: run_wann()
elif args.crpa: run_crpa()
else: run_vasp()
