#!/usr/bin/env python3

"""
Project: DFTHub - First-Principles Hubbard Parameter Extraction
Description:
	Extracts Hubbard model parameters (t, U, J) based on first-principles
	calculations performed with ABINIT.
	- Hopping parameters (t) are derived from Wannier90 tight-binding
	  Hamiltonian constructed from ABINIT bloch states.
	- Interaction parameters (U, J) are extracted from constrained Random
	  Phase Approximation (cRPA) calculations within the ABINIT framework.

Author: Yerin Jang
Affiliation: Gwangju Institute of Science and Technology, Department of Materials Science and Engineering
Last Modified: 2026 . 04 . 03
"""

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

def clean(target_dir, target_type):
	if not target_type in ['cell', 'band', 'wann', 'crpa']:
		print(f'\'{target_type}\' is wrong type')
		sys.exit(1)

	for fn in [f'{target_dir}/{fn}' for fn in os.listdir(target_dir) if re.search(f'{target_type}_', fn)]:
		if not re.search('[.]abi|[.]win', fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

def get_density(target_dir, band_idx=2, is_sppol=False):
	target_mat = target_dir.split('-')[0]

	if is_sppol:
		print('sppol is not allowed')
		sys.exit(1)

	# ABINIT: [root]_DEN.cube
	fn = f'{target_dir}/{target_mat}_wann_o_DS1_DEN'

	### convert ABINIT binary density to Gaussian Cube format
	den2cube = subprocess.run(
			['cut3d'],
			input = '\n'.join([fn, '14', fn+'.cube', '0', '0']),
			text=True,
			capture_output=True
			)
	if den2cube.returncode != 0:
		print(den2cube.stderr)
		sys.exit(1)

	with open(fn+'.cube', 'r') as f: lines = f.readlines()
	ls = [lines[i].split() for i in [2, 3, 4, 5]]
	natoms = int(ls[0][0])
	O = np.array(ls[0][1:], dtype=float) # Origin
	N = np.array((ls[1][0], ls[2][0], ls[3][0]), dtype=int) # Real-space grid number
	#S = np.array((ls[1][1:], ls[2][1:], ls[3][1:]), dtype=float) # Real-space grid vector
	s = np.meshgrid(np.arange(N[0])/N[0], np.arange(N[1])/N[1], np.arange(N[2])/N[2], indexing='ij') # Fractional real-space grid vector
	rho_dft = np.array(lines[6+natoms:], dtype=float).reshape(N) # DFT density in real space (rho_DFT(r))

	# ABINIT: [root]_GSR.nc
	fn = f'{target_dir}/{target_mat}_wann_o_DS1_GSR.nc'
	with abilab.abiopen(abidata.ref_file(f'{os.getcwd()}/{fn}')) as f:
		fk = f.ebands.occfacts[0, :, 1]
		wk = f.ebands.kpoints.weights	
	
	# WANNIER90: [root]_u.mat
	fn = f'{target_dir}/w90_u.mat'
	with open(fn, 'r') as f: lines = f.readlines()
	nkpts, nwann, _ = map(int, lines[1].split())
	kpts = np.loadtxt(lines[3::3])
	U = np.loadtxt(lines[4::3])
	U = (U[:, 0] + 1j*U[:, 1]).reshape(nkpts, nwann, nwann)

	# WANNIER90: [root]_u_dis.mat
	# If disentanglement exist, final U shape is (nkpts, nband, nwann)

	# WANNIER90: [root]_hr.dat
	fn = f'{target_dir}/w90_hr.dat'
	with open(fn, 'r') as f: lines = f.readlines()
	nrpts = int(lines[2])
	rpts = np.unique(np.array([l.split()[:3] for l in lines[-nrpts*nwann*nwann:]], dtype=int), axis=0)
	rpts_dict = {tuple(r): i for i, r in enumerate(rpts)}

	# ABINIT: UNKp.s
	phi = np.zeros((nwann, nrpts, *N), dtype=complex) # Wannier function in real space (phi_iR(r))
	for k in range(nkpts):
		with scipy.io.FortranFile(f'{target_dir}/UNK{k+1:05d}.1', 'r') as f:
			f.read_record(np.int32)
			u = f.read_record(np.complex128).reshape(N, order='F')
		phase_s = np.exp(1j * 2*np.pi * (kpts[k, 0]*s[0] + kpts[k, 1]*s[1] + kpts[k, 2]*s[2]))
		psi = phase_s * u
		for r in range(nrpts):
			for i in range(nwann):
				phase_r = np.exp(-1j * 2*np.pi * np.dot(kpts[k], rpts[r]))
				phi[i, r] += wk[k] * phase_r * U[k, 0, i] * psi

	"""
	### Wannier function visualization
	# phase alignment
	r0_idx = np.where((rpts == [0, 0, 0]).all(axis=1))[0][0]
	max_idx = np.unravel_index(np.argmax(np.abs(phi[0, r0_idx])), phi[0, r0_idx].shape)
	theta = np.angle(phi[0, r0_idx][max_idx])
	phi *= np.exp(-1j * theta) 

	phi = np.roll(phi, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(2, 3, 4))
	plt.imshow(np.real(phi[0, r0_idx][N[0]//2, :, :]), origin='lower')
	plt.colorbar()
	plt.show()

	# WANNIER90: [root].cube 
	fn = f'{target_dir}/w90_00001.cube'
	with open(fn, 'r') as f: lines = f.readlines()
	phi = np.array([x for line in lines[7:] for x in line.split()], dtype=float).reshape(94, 95, 95)
	plt.imshow(np.real(phi[94//2, :, :]), origin='lower')
	plt.colorbar()
	plt.show()
	"""

	# Occupancy matrix from DFT in Wannier basis (<c^dagger_iR c_jR'>_DFT)
	n_dft = np.zeros((nwann, nwann, nrpts), dtype=complex)
	for k in range(nkpts):
		#occ = U[k].conj().T @ fk[k] @ U[k]
		occ =  fk[k] * np.outer(U[k, 0, :].conj().T, U[k, 0, :])
		for r in range(nrpts):
			phase_r = np.exp(1j * 2*np.pi * np.dot(kpts[k], rpts[r]))
			n_dft[:, :, r] += wk[k] * phase_r * occ

	# Occupancy matrix from QC (<c^dagger_iR c_jR'>_QC)
	n_qc = 0.9 * np.ones((nwann, nwann, nrpts), dtype=complex)

	# Corrected occupancy matrix (<c^dagger_iR c_jR'>_QC - <c^dagger_iR c_jR'>_DFT)
	n_corr = (n_qc - n_dft).real

	rho_corr = np.zeros(np.prod(N), dtype=complex)
	for r1 in range(nrpts):
		phi1 = phi[:, r1].reshape(nwann, -1)
		for r2 in range(nrpts):
			phi2 = phi[:, r2].reshape(nwann, -1)
			dr_idx = rpts_dict.get(tuple(rpts[r2] - rpts[r1]))
			if not dr_idx is None:
				rho_corr += np.sum(phi1.conj() * (n_corr[:, :, dr_idx] @ phi2), axis=0)
	rho_corr = rho_corr.reshape(N).real

	rho = rho_dft + rho_corr

	"""
	### total density visualization
	rho = np.roll(rho, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(0, 1, 2))
	plt.imshow(rho[N[0]//2, :, :], origin='lower')
	plt.colorbar()
	plt.show()
	"""

if args.clean: clean(args.target_dir, args.clean)
elif args.density: get_density(args.target_dir, band_idx=2, is_sppol=False)
else:
	parser.print_help()
	sys.exit(1)
