#!/usr/bin/env python3

import os
import re
import sys
import time
import shutil
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
plt.rcParams['font.size'] = 20

from ase import Atoms, io
from ase.dft.kpoints import bandpath, get_special_points

from abipy import abilab
import abipy.data as abidata

parser = argparse.ArgumentParser()
parser.add_argument('mat', type=str)
parser.add_argument('--clean', type=str)
parser.add_argument('--poscar', action='store_true')
parser.add_argument('--hist', action='store_true')
parser.add_argument('--etot', action='store_true')
parser.add_argument('--kpts', action='store_true')
parser.add_argument('--wann', action='store_true')
parser.add_argument('--band', type=int, nargs='?', const=-1, default=99)
args = parser.parse_args()

matdir = f'{os.getcwd()}/{args.mat}'
mat = args.mat.split('-')[0]

BOHR_TO_ANG = 0.52917721092
HA_TO_EV = 27.211386245

def clean(target):
	for fn in os.listdir(matdir):
		fn = f'{matdir}/{fn}'
		if re.search(f'{target}_o|{target}_.abo', fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

def gen_poscar(cell=[0, 0, 0], size=(1, 1, 1), fn='POSCAR'):
	fn = f'{matdir}/{fn}'
	bulk = Atoms(mat, cell=cell, pbc=True) * size
	io.write(fn, bulk, 'vasp')

	print(f'{fn} ({bulk.cell.get_bravais_lattice()}):')
	with open(fn, 'r') as f: print(f.read())
	print('* Please add element symbol after the atomic position')

def gen_hspts(path, nkpts):
	fn = f'{matdir}/POSCAR'
	bulk = io.read(fn)

	kpts, hspts, hspts_label = bandpath(path, bulk.cell, npoints=nkpts).get_linear_kpoint_axis()
	hspts = [np.where(np.isclose(kpts, p, atol=1e-6))[0][0] for p in hspts]
	hspts_label = [l if l != 'G' else r'$\Gamma$' for l in hspts_label]

	return bulk, hspts, hspts_label

def show_hist():
	fn = f'{matdir}/{mat}_cell_o_HIST.nc'
	with abilab.abiopen(abidata.ref_file(fn)) as f:
		print(f'{fn}:')
		for i, (struct, etot) in enumerate(zip(f.structures, f.etotals)):
			print(f'itr: {i}  Volume: {struct.volume:.6f}  Etotal: {etot:.6f}')
		print('\n', f.final_structure, sep='', end='\n\n')

		#cell = np.round(f.final_structure.lattice.abc, decimals=6)
		#gen_poscar(cell=cell, size=(1, 1, 1), fn='POSCAR')
		#gen_poscar(cell=cell, size=(2, 2, 1), fn='POSCAR2')
	
def show_etot():
	fin, acell, etot = 0, [], []
	with open(f'{matdir}/{mat}_cell_.abo', 'r') as f:
		for line in f:
			if fin:
				if 'acell' in line: acell.append(line.split()[1])
				if 'etotal' in line: etot.append(line.split()[1])
			elif 'END DATASET' in line: fin = 1
	data = np.column_stack((acell[:-1], etot)).astype(float)
	data = np.column_stack((data[:, 0] * BOHR_TO_ANG, data[:, 1] * HA_TO_EV))

	fig, ax = plt.subplots(figsize=(8, 6), dpi=120, tight_layout=True)
	ax.plot(data[:, 0], data[:, 1], '.-')
	ax.axvline(data[np.argmin(data[:, 1]), 0], color='r', ls='--')
	ax.set_xticks(data[:, 0], [f'{x:.2f}' for x in data[:, 0]], rotation=60, ha='center')
	ax.set_xlabel(r'$d$ (Å)')
	ax.set_ylabel(r'$E_{\mathrm{tot}}$ (eV)')

	fig.savefig(f'{matdir}/fig/{mat}_cell_etot.svg')
	plt.show()
	
def show_kpts(path, nkpts):
	bulk, hspts, hspts_label = gen_hspts(path, nkpts)

	ndivk = [hspts[i] - hspts[i-1] for i in range(1, len(hspts))]
	print('ndivk:\n', *ndivk)

	kptbounds = get_special_points(bulk.cell)
	print('\nkptbounds:')
	for kpt in path:
		print(*kptbounds[kpt], f'#{kpt}')

def show_fatband(l, nkpts, nband, hspts, hspts_label, ylim=(-9, 9)):
	fn_head = f'{mat}_band3_o_DS2_FATBANDS_at0001_{mat}_is1'
	#fn_head = f'{mat}_band_o_DS2_FATBANDS_at0001_{mat}'
	if l < 0:
		fn_list = sorted([fn for fn in os.listdir(matdir) if re.search(fn_head, fn)])
	else:
		fn_list = [f'{fn_head}_l{l:04d}']
	print(*fn_list, sep='\n', end='\n\n')

	band = []
	for fn in fn_list:
		with open(f'{matdir}/{fn}', 'r') as f:
			band_l = []
			for line in f:
				if line.strip() and not re.match(r'[#@&]', line):
					band_l.append([float(v) for v in line.strip().split()])
			band.append(np.split(np.array(band_l), nband))

	fig, ax = plt.subplots(figsize=(8, 6), dpi=120, tight_layout=True)
	ecolor = list(mcolors.TABLEAU_COLORS)

	L = len(band) // 2
	w, alpha = 6, 0.7
	for l in range(2):
		for i, (b1, b2) in enumerate(zip(band[l], band[l+L])):
			ax.errorbar(b1[:, 0], b1[:, 1], b1[:, 2]*w, c='k', ecolor=ecolor[l], alpha=alpha, label=None if i else rf'$l = {l}$')
			ax.errorbar(b2[:, 0], b2[:, 1], b2[:, 2]*w, c='k', ecolor=ecolor[l], alpha=alpha)
	ax.axhline(0, lw=2, ls='--', color='k', alpha=0.2)

	for p in hspts:
		ax.axvline(p, lw=2, color='k', alpha=0.2)
	ax.set_xticks(hspts, labels=hspts_label)

	ax.set_xlim(0, nkpts-1)
	ax.set_ylim(*ylim)
	ax.set_ylabel(r'$E - E_{F}$')
	ax.legend(fontsize='x-small', loc='upper right')

	#dos = []
	#with open(f'{matdir}/{mat}_band_o_DS3_DOS_AT0001', 'r') as f:
	#	for line in f:
	#		if 'Fermi energy' in line: fermi = float(re.search(r'[-]?\d+\.\d+', line).group())
	#		if line.strip() and not re.match(r'[#@&]', line):
	#			dos.append([float(v) for v in line.strip().split()])
	#dos = np.split(np.array(dos), 2)

	#fig, ax = plt.subplots(1, 2, figsize=(9, 6), width_ratios=(3, 1), dpi=120, tight_layout=True)
	#ecolor = list(mcolors.TABLEAU_COLORS)

	#L = len(band) // 2
	#w, alpha = 6, 0.7
	#for l in range(2):
	#	for i, (b1, b2) in enumerate(zip(band[l], band[l+L])):
	#		ax[0].errorbar(b1[:, 0], b1[:, 1], b1[:, 2]*w, c='k', ecolor=ecolor[l], alpha=alpha, label=None if i else rf'$l = {l}$')
	#		ax[0].errorbar(b2[:, 0], b2[:, 1], b2[:, 2]*w, c='k', ecolor=ecolor[l], alpha=alpha)
	#		#if b1[0, 1] > ylim[0] and b1[0, 1] < ylim[1]: ax[0].text(b1[0, 0]+1, b1[0, 1]+0.1, rf'{i+1}$\uparrow$', fontsize='xx-small')
	#		#if b2[0, 1] > ylim[0] and b2[0, 1] < ylim[1]: ax[0].text(b2[0, 0]+1, b2[0, 1]+0.1, rf'{i+1}$\downarrow$', fontsize='xx-small')
	#ax[0].axhline(0, lw=2, ls='--', color='k', alpha=0.2)

	#for dos_s, label, color in zip(dos, [r'$s\uparrow$', r'$s\downarrow$'], ['r', 'b']):
	#	ax[1].plot(dos_s[:, 1], (dos_s[:, 0] - fermi)*HA_TO_EV, label=label, c=color)
	#ax[1].axhline(0, lw=2, ls='--', color='k', alpha=0.2)

	##_, hspts, hspts_label = gen_hspts(path, nkpts)
	##for p in hspts[1:-1]:
	#for p in hspts:
	#	ax[0].axvline(p, lw=2, color='k', alpha=0.2)
	#ax[0].set_xticks(hspts, labels=hspts_label)

	#ax[0].set_xlim(0, nkpts-1)
	#ax[0].set_ylim(*ylim)
	#ax[0].set_ylabel(r'$E - E_{F}$')
	#ax[0].legend(fontsize='x-small', loc='upper right')

	#ax[1].set_ylim(*ylim)
	#ax[1].legend(fontsize='x-small', loc='upper right')

	fig.savefig(f'{matdir}/fig/{fn_head if l < 0 else fn_list[0]}.svg')
	plt.show()

def show_wann(nkpts, nband, hspts, hspts_label, ylim):
	fn_abo = f'{mat}2_wann_.abo'
	#fn_wan1 = f'w90_up_band.dat'
	#fn_wan2 = f'w90_down_band.dat'
	fn_wan = f'w90_band.dat'
	#fn_dft1 = f'{mat}_band_o_DS2_FATBANDS_at0001_{mat}_is1_l0000'
	#fn_dft2 = f'{mat}_band_o_DS2_FATBANDS_at0001_{mat}_is2_l0000'
	fn_dft = f'{mat}2_band_o_DS2_FATBANDS_at0001_{mat}_is1_l0000'

	with open(f'{matdir}/{fn_abo}', 'r') as f:
		for line in f:
			if 'Fermi' in line: fermi = float(re.search(r'[-]?\d+\.\d+', line).group()) * HA_TO_EV

	band_wan = []
	#for fn_wan in [fn_wan1, fn_wan2]:
	for fn_wan in [fn_wan]:
		with open(f'{matdir}/{fn_wan}', 'r') as f:
			band_wan_s = []
			for line in f:
				if line.strip():
					band_wan_s.append([float(v) for v in line.strip().split()[:2]])
		band_wan.append(np.array(band_wan_s))

	band_dft = []
	#for fn_dft in [fn_dft1, fn_dft2]:
	for fn_dft in [fn_dft]:
		with open(f'{matdir}/{fn_dft}', 'r') as f:
			band_dft_s = []
			for line in f:
				if line.strip() and not re.match(r'[#@&]', line):
					band_dft_s.append([float(v) for v in line.strip().split()])
			band_dft.append(np.split(np.array(band_dft_s), nband))

	fig, ax = plt.subplots(figsize=(8, 6), dpi=120, tight_layout=True)
	#for s in range(2):
	for s in range(1):
		ax.scatter(np.arange(len(band_wan[s])), band_wan[s][:, 1]-fermi, color='r', label=None if s else 'wann')
		for b in band_dft[s]:
			ax.plot(b[:, 0], b[:, 1], c='k')
	ax.axhline(0, lw=2, ls='--', color='k', alpha=0.2)

	#_, hspts, hspts_label = gen_hspts(path, nkpts)
	#for p in hspts[1:-1]:
	for p in hspts:
		ax.axvline(p, lw=2, color='k', alpha=0.2)
	ax.set_xticks(hspts, labels=hspts_label)

	ax.set_xlim(0, nkpts-1)
	ax.set_ylim(*ylim)
	ax.set_ylabel(r'$E - E_{F}$')
	ax.legend(fontsize='x-small', loc='upper right')

	fig.savefig(f'{matdir}/fig/wan_band.svg')
	plt.show()
		
if args.clean: clean(args.clean)
elif args.poscar:
	if   re.match('Li', mat): gen_poscar(cell=[20., 20., 20.])
	elif re.match('Pt', mat): gen_poscar(cell=[2.79, 2.79, 20.])
elif args.hist: show_hist()
elif args.etot: show_etot()
elif args.kpts: show_kpts(path='GX', nkpts=100)
elif args.wann:
	if re.match('Li', mat): show_wann(nkpts=100, nband=8, hspts=[0, 99], hspts_label=[r'$\Gamma$', 'X'], ylim=(-1.5, 3.5))
elif args.band < 99:
	if re.match('Li', mat): show_fatband(args.band, nkpts=100, nband=8, hspts=[0, 99], hspts_label=[r'$\Gamma$', 'X'], ylim=(-1.5, 3.5))
else:
	parser.print_help()
	sys.exit(1)
