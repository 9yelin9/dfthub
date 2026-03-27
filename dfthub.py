#!/usr/bin/env python3

import os
import re
import sys
import time
import scipy
import shutil
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt

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

def get_density(target_dir, band_idx=2):
	target_mat = target_dir.split('-')[0]

	if is_sppol:
		print('sppol is not allowed')
		sys.exit(1)

	# ABINIT: *_DEN.cube
	fn = f'{target_dir}/{target_mat}_wann_o_DS1_DEN'
	"""
	den2cube = subprocess.run(
			['cut3d'],
			input = '\n'.join([fn, '14', fn+'.cube', '0', '0']),
			text=True,
			capture_output=True
			)
	if den2cube.returncode != 0:
		print(den2cube.stderr)
		sys.exit(1)
	"""

	with open(fn+'.cube', 'r') as f: lines = f.readlines()
	ls = [lines[i].split() for i in [2, 3, 4, 5]]
	natoms = int(ls[0][0])
	O = np.array(ls[0][1:], dtype=float)
	N = np.array((ls[1][0], ls[2][0], ls[3][0]), dtype=int)
	A = np.array((ls[1][1:], ls[2][1:], ls[3][1:]), dtype=float)
	rho_dft = np.reshape(lines[6+natoms:], N).astype(float)

	# ABINIT: *_GSR.nc
	fn = f'{target_dir}/{target_mat}_wann_o_DS1_GSR.nc'
	with abilab.abiopen(abidata.ref_file(f'{os.getcwd()}/{fn}')) as f:
		occ = f.ebands.occfacts[0, :, 1]
		w = f.ebands.kpoints.weights	
	
	# WANNIER90: *_u.mat
	fn = f'{target_dir}/w90_u.mat'
	with open(fn, 'r') as f: lines = f.readlines()
	nkpts, nwann, _ = map(int, lines[1].split())
	kpts = np.loadtxt(lines[3::3])
	U = np.loadtxt(lines[4::3])
	U = np.reshape(U[:, 0] + 1j*U[:, 1], (nkpts, nwann, nwann))

	# ABINIT: UNK*
	s = np.meshgrid(np.arange(N[0])/N[0], np.arange(N[1])/N[1], np.arange(N[2])/N[2], indexing='ij')
	phi = np.zeros(N, dtype=complex)
	for k in range(nkpts):
		with scipy.io.FortranFile(f'{target_dir}/UNK{k+1:05d}.1', 'r') as f:
			f.read_record(np.int32)
			u = np.reshape(f.read_record(np.complex128), N, order='F')
		phase = np.exp(1j * 2*np.pi * (kpts[k, 0]*s[0] + kpts[k, 1]*s[1] + kpts[k, 2]*s[2]))
		psi = phase * u
		phi += w[k] * U[k, 0, 0] * psi
	max_idx = np.unravel_index(np.argmax(np.abs(phi)), phi.shape)
	theta = np.angle(phi[max_idx])
	phi *= np.exp(-1j * theta) 

	"""
	phi = np.roll(phi, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(0, 1, 2))
	plt.imshow(np.real(phi[N[0]//2, :, :]), origin='lower')
	plt.colorbar()
	plt.show()

	# WANNIER90: *.cube
	fn = f'{target_dir}/w90_00001.cube'
	with open(fn, 'r') as f: lines = f.readlines()
	phi = np.reshape([x for line in lines[7:] for x in line.split()], (94, 95, 95)).astype(float)
	plt.imshow(np.real(phi[94//2, :, :]), origin='lower')
	plt.colorbar()
	plt.show()
	"""

	n_qc = 0.9 * np.ones((nwann, nwann), dtype=complex)
	n_dft = np.zeros((nwann, nwann), dtype=complex)
	for k in range(nkpts):
		n_dft += w[k] * (U[k, 0, 0].conj() * occ[k] * U[k, 0, 0])
	n_qc, n_dft = n_qc.real, n_dft.real

	rho_corr = (np.abs(phi)**2) * (n_qc - n_dft)
	rho = rho_dft + rho_corr

	rho = np.roll(rho, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(0, 1, 2))
	plt.imshow(rho[N[0]//2, :, :], origin='lower')
	plt.colorbar()
	plt.show()

if args.clean: clean(args.target_dir, args.clean)
elif args.density: get_density(args.target_dir)
else:
	parser.print_help()
	sys.exit(1)
