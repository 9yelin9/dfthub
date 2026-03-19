#!/usr/bin/env python3

import os
import re
import sys
import time
import shutil
import argparse
import subprocess
import numpy as np

from abipy import abilab
import abipy.data as abidata

parser = argparse.ArgumentParser()
parser.add_argument('target_dir', type=str)
parser.add_argument('-c', '--clean', type=str)
parser.add_argument('-d', '--density', action='store_true')
args = parser.parse_args()

is_sppol = True if args.target_dir in ['Li'] else False

def clean(target_dir, target_type):
	if not target_type in ['cell', 'band', 'wann', 'crpa']:
		print(f'\'{target_type}\' is wrong type')
		sys.exit(1)

	for fn in [f'{target_dir}/{fn}' for fn in os.listdir(target_dir) if re.search(f'{target_type}_', fn)]:
		if not re.search('[.]abi|[.]win', fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

def read_cube(lines):
	l2 = lines[2].split()
	l3 = lines[3].split()
	l4 = lines[4].split()
	l5 = lines[5].split()

	cube = {
			'natoms': int(l2[0]),
			'origin': np.array(l2[1:], dtype=float),
			'nx': int(l3[0]),
			'ny': int(l4[0]),
			'nz': int(l5[0]),
			'ax': np.array(l3[1:], dtype=float),
			'ay': np.array(l4[1:], dtype=float),
			'az': np.array(l5[1:], dtype=float),
			}
	return cube

def get_density(target_dir):
	target_mat = target_dir.split('-')[0]

	if is_sppol:
		print('sppol is not allowed')
		sys.exit(1)

	fd = f'{target_dir}/{target_mat}_wann_o_DS1_DEN.cube'
	fw = f'{target_dir}/w90_00001.cube'
	fu = f'{target_dir}/w90_u.mat'
	ff = f'{target_dir}/{target_mat}_wann_o_DS1_GSR.nc'

	"""
	den2cube = subprocess.run(
			['cut3d'],
			input = '\n'.join([re.sub('[.]cube', '', fd), '14', fd, '0', '0']),
			text=True,
			capture_output=True
			)
	if den2cube.returncode != 0:
		print(den2cube.stderr)
		sys.exit(1)
	"""

	"""
	with open(fd, 'r') as f: lines = f.readlines()
	dcube = read_cube(lines)
	rho = np.reshape(lines[6+dcube['natoms']:], (dcube['nx'], dcube['ny'], dcube['nz'])).astype(float)
	
	with open(fw, 'r') as f: lines = f.readlines()
	wcube = read_cube(lines)
	phi = np.reshape([x for line in lines[6+wcube['natoms']:] for x in line.split()],
			(wcube['nx'], wcube['ny'], wcube['nz'])).astype(float)

	"""

	with open(fu, 'r') as f: lines = [lines.strip() for lines in f.readlines() if lines.strip()]
	nkpts, nwann, _ = map(int, lines[1].split())
	kpts = np.loadtxt(lines[2::2])
	U = np.loadtxt(lines[3::2])
	U = np.reshape(U[:, 0] + 1j*U[:, 1], (nkpts, nwann, nwann))

	with abilab.abiopen(abidata.ref_file(f'{os.getcwd()}/{ff}')) as f:
		occ = f.ebands.occfacts[0, :, 1]
		wk = f.ebands.kpoints.weights	

	n = np.zeros((nwann, nwann), dtype=complex)
	for k in range(nkpts):
		#n += wk[k] * (U[k].conj().T @ occ[k] @ U[k])
		n += wk[k] * occ[k]
	print(n)

if args.clean: clean(args.target_dir, args.clean)
elif args.density: get_density(args.target_dir)
else:
	parser.print_help()
	sys.exit(1)
