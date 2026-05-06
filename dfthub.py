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
import h5py
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
parser.add_argument('-d', '--density', type=str)
args = parser.parse_args()

def resub_float(string_sub, string_search):
	return float(re.sub(string_sub, '', re.search(rf'{string_sub}[-]?\d+[.]\d+', string_search).group()))

def clean(target_dir, target_type):
	if not target_type in ['cell', 'band', 'wann', 'crpa']:
		print(f'\'{target_type}\' is wrong type')
		sys.exit(1)

	for fn in [f'{target_dir}/{fn}' for fn in os.listdir(target_dir) if re.search(f'{target_type}_', fn)]:
		if not re.search('[.]abi|[.]win', fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

def get_density(target_dir, cdc_qc_file, nn=2, show_digit=3, is_sppol=False):
	print(f'\ntarget_dir={target_dir}')
	print(f'cdc_qc_file={cdc_qc_file}')
	print(f'nn={nn}')
	print(f'is_sppol={is_sppol}\n')

	if is_sppol:
		print('sppol is not allowed')
		sys.exit(1)

	root_tag = args.target_dir.split('_')[0]
	base_dir = f'{root_tag}_base'
	save_dir = f'{root_tag}_t{resub_float('t', cdc_qc_file):.3f}_tp{resub_float('tp', cdc_qc_file):.3f}_U{resub_float('U', cdc_qc_file):.3f}'
	os.makedirs(save_dir, exist_ok=True)
	for fn in ['wann_.abi', 'crpa_.abi', 'w90.win']:
		src = os.path.join(base_dir, fn)
		dst = os.path.join(save_dir, fn)
		shutil.copy(src, dst)
	print(f'save_dir={save_dir}\n')

	# ABINIT: [root]_DEN.cube
	fn = f'{target_dir}/wann_o_DS1_DEN.cube'
	if not os.path.isfile(fn):
		den2cube = subprocess.run( # convert ABINIT binary density to Gaussian Cube format
				['cut3d'],
				input = '\n'.join([re.sub('.cube', '', fn), '14', fn, '0', '0']),
				text=True,
				capture_output=True
				)
		if den2cube.returncode != 0:
			print(den2cube.stderr)
			sys.exit(1)

	with open(fn, 'r') as f: lines = f.readlines()
	data = [lines[i].split() for i in [2, 3, 4, 5]]
	natoms = int(data[0][0])
	O = np.array(data[0][1:], dtype=float) # Origin
	N = np.array((data[1][0], data[2][0], data[3][0]), dtype=int) # Real-space grid number
	S = np.array((data[1][1:], data[2][1:], data[3][1:]), dtype=float) # Real-space grid vector
	dV = np.abs(np.dot(S[0], np.cross(S[1], S[2])))
	s = np.meshgrid(np.arange(N[0])/N[0], np.arange(N[1])/N[1], np.arange(N[2])/N[2], indexing='ij') # Fractional real-space grid vector
	rho_dft = np.array(lines[6+natoms:], dtype=float).reshape(N) # DFT density in real space (rho_DFT(r))
	print(f'{fn}:')
	print(f'grid num={N}')
	print(f'rho_dft sum={np.round(np.sum(rho_dft) * dV, show_digit)}\n')

	# ABINIT: [root]_GSR.nc
	fn = f'{target_dir}/wann_o_DS1_GSR.nc'
	with abilab.abiopen(abidata.ref_file(f'{os.getcwd()}/{fn}')) as f:
		fk = f.ebands.occfacts[0, :, 1]
		wk = f.ebands.kpoints.weights	
	print(f'{fn}:')
	print(f'fk={np.round(fk, show_digit)}')
	print(f'wk={np.round(wk, show_digit)}\n')
	
	# WANNIER90: [root]_u.mat
	fn = f'{target_dir}/w90_u.mat'
	with open(fn, 'r') as f: lines = f.readlines()
	nkpts, nwann, _ = map(int, lines[1].split())
	kpts = np.loadtxt(lines[3::3])
	U = np.loadtxt(lines[4::3])
	U = (U[:, 0] + 1j*U[:, 1]).reshape(nkpts, nwann, nwann)
	print(f'{fn}:')
	print(f'nkpts={nkpts}')
	print(f'nwann={nwann}')
	print(f'U norm={np.round([np.linalg.norm(U[k]) for k in range(nkpts)], show_digit)}\n')

	# WANNIER90: [root]_u_dis.mat
	# If disentanglement exist, final U shape is (nkpts, nband, nwann)

	# WANNIER90: [root]_hr.dat
	fn = f'{target_dir}/w90_hr.dat'
	with open(fn, 'r') as f: lines = f.readlines()
	nrpts0 = int(lines[2])
	rpts0 = np.array([l.split()[:3] for l in lines[-nrpts0*nwann*nwann:]], dtype=int)
	dist = np.linalg.norm(rpts0, axis=1)
	rpts = rpts0[np.isin(dist, np.sort(np.unique(dist))[:nn+1])]
	rpts_dict = {tuple(r): i for i, r in enumerate(rpts)}
	nrpts = len(rpts)
	print(f'{fn}:')
	print(f'nrpts={nrpts}')
	print(f'rpts=\n{rpts}\n')

	# ABINIT: UNKp.s
	fn = f'{target_dir}/UNK'
	dv = 1 / np.sqrt(np.prod(N))
	phi = np.zeros((nwann, nrpts, *N), dtype=complex) # Wannier function in real space (phi_iR(r))
	for k in range(nkpts):
		with scipy.io.FortranFile(f'{fn}{k+1:05d}.1', 'r') as f:
			f.read_record(np.int32)
			u = f.read_record(np.complex128).reshape(N, order='F') * dv
		phase_s = np.exp(1j * 2*np.pi * (kpts[k, 0]*s[0] + kpts[k, 1]*s[1] + kpts[k, 2]*s[2]))
		psi = phase_s * u
		for r in range(nrpts):
			for i in range(nwann):
				phase_r = np.exp(-1j * 2*np.pi * np.dot(kpts[k], rpts[r]))
				phi[i, r] += wk[k] * phase_r * U[k, 0, i] * psi
	phi_norm = [np.sum(np.abs(phi[0, r])**2) for r in range(nrpts)]
	print(f'{fn}:')
	print(f'|phi[i, r]|^2={np.round(phi_norm, show_digit)} -> sum={np.round(np.sum(phi_norm), show_digit)}\n')

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
	cdc_dft = np.zeros((nwann, nwann, nrpts), dtype=complex)
	for k in range(nkpts):
		#occ = U[k].conj().T @ fk[k] @ U[k]
		occ =  fk[k] * np.outer(U[k, 0, :].conj().T, U[k, 0, :])
		for r in range(nrpts):
			phase_r = np.exp(1j * 2*np.pi * np.dot(kpts[k], rpts[r]))
			cdc_dft[:, :, r] += wk[k] * phase_r * occ
	cdc_dft = cdc_dft.real

	# Occupancy matrix from QC (<c^dagger_iR c_jR'>_QC)
	with h5py.File(cdc_qc_file, 'r') as f:
		cdc_up = f['cdagc_up'][:]
		cdc_dn = f['cdagc_dn'][:]
	mid = cdc_up.shape[0] // 2
	cdc_qc = (cdc_up[mid, mid-nn:mid+nn+1] + cdc_dn[mid, mid-nn:mid+nn+1]).reshape(nwann, nwann, nrpts)

	# Corrected occupancy matrix (<c^dagger_iR c_jR'>_QC - <c^dagger_iR c_jR'>_DFT)
	cdc_cor = cdc_qc - cdc_dft
	print(f'cdc_qc ={np.round(cdc_qc,  show_digit)}')
	print(f'cdc_dft={np.round(cdc_dft, show_digit)}')
	print(f'cdc_cor={np.round(cdc_cor, show_digit)}\n')

	rho_cor = np.zeros(np.prod(N), dtype=complex)
	for r1 in range(nrpts):
		phi1 = phi[:, r1].reshape(nwann, -1)
		for r2 in range(nrpts):
			phi2 = phi[:, r2].reshape(nwann, -1)
			dr_idx = rpts_dict.get(tuple(rpts[r2] - rpts[r1]))
			if not dr_idx is None:
				rho_cor += np.sum(phi1.conj() * (cdc_cor[:, :, dr_idx] @ phi2), axis=0)
	rho_cor = rho_cor.reshape(N).real
	rho_new = rho_dft + rho_cor
	print(f'rho_dft={np.round([rho_dft.min(), rho_dft.max()], show_digit)} -> sum={np.round(np.sum(rho_dft) * dV, show_digit)}')
	print(f'rho_cor={np.round([rho_cor.min(), rho_cor.max()], show_digit)} -> sum={np.round(np.sum(rho_cor) * dV, show_digit)}')
	print(f'rho_new={np.round([rho_new.min(), rho_new.max()], show_digit)} -> sum={np.round(np.sum(rho_new) * dV, show_digit)}\n')

	# Save new density
	rho_new_bytes = rho_new.flatten(order='F').astype(np.float64).tobytes()
	fn_dft = f'{target_dir}/wann_o_DS1_DEN'
	fn_new = f'{save_dir}/dfthub_DEN'
	shutil.copyfile(fn_dft, fn_new)
	with open(fn_new, 'r+b') as f:
		f.seek(0, 2)
		fsize = f.tell()
		f.seek(fsize - len(rho_new_bytes) - 4, 2)
		f.write(rho_new_bytes)
	print(f'{fn_new} generated\n')

	"""
	### Density visualization
	rho_dft_roll = np.roll(rho_dft, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(0, 1, 2))[N[0]//2, :, :]
	rho_new_roll = np.roll(rho_new, shift=(N[0]//2, N[1]//2, N[2]//2), axis=(0, 1, 2))[N[0]//2, :, :]
	vmin, vmax = rho_new_roll.min(), rho_new_roll.max()

	fig, ax = plt.subplots(1, 2, figsize=(8, 6), constrained_layout=True)
	im0 = ax[0].imshow(rho_dft_roll, origin='lower', vmin=vmin, vmax=vmax)
	im1 = ax[1].imshow(rho_new_roll, origin='lower', vmin=vmin, vmax=vmax)
	ax[0].set_title('rho_dft')
	ax[1].set_title('rho_new')
	fig.colorbar(im1, ax=ax, shrink=0.5)
	plt.show()
	"""

if args.clean: clean(args.target_dir, args.clean)
elif args.density: get_density(args.target_dir, args.density)
else:
	parser.print_help()
	sys.exit(1)
