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

parser = argparse.ArgumentParser()
parser.add_argument('mat', type=str)
<<<<<<< HEAD
<<<<<<< HEAD
parser.add_argument('--clean', type=str, choices=['cell', 'band', 'wann', 'crpa'])
=======
parser.add_argument('--clean', action='store_true')
>>>>>>> parent of 1d84223 (Update)
parser.add_argument('--poscar', action='store_true')
parser.add_argument('--kpts', action='store_true')
<<<<<<< HEAD
parser.add_argument('--wann', action='store_true')
parser.add_argument('--band', type=int, nargs='?', const=-1, default=99)
=======
parser.add_argument('--band', type=int, nargs='+', default=-1)
>>>>>>> parent of 1d84223 (Update)
=======
parser.add_argument('--clean', action='store_true')
parser.add_argument('--poscar', action='store_true')
parser.add_argument('--kpts', action='store_true')
parser.add_argument('--band', type=int, nargs='+', default=-1)
>>>>>>> parent of 1d84223 (Update)
args = parser.parse_args()

matdir = f'{os.getcwd()}/{args.mat}'
poscar = f'{matdir}/POSCAR'

def clean(keep=r'POSCAR|.*\.abi|fig'):
	for fn in os.listdir(matdir):
		fn = f'{matdir}/{fn}'
		if not re.search(keep, fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

<<<<<<< HEAD
<<<<<<< HEAD
def gen_poscar(fn, cell, size):
	bulk = Atoms('Pt', cell=cell, pbc=True) * size
	
	fn = f'{matdir}/{fn}'
	io.write(fn, bulk, 'vasp')
=======
def gen_poscar(d):
	bulk = Atoms('Pt', cell=[d, d, 10.], pbc=True)
	io.write(poscar, bulk, 'vasp')
>>>>>>> parent of 1d84223 (Update)
=======
def gen_poscar(d):
	bulk = Atoms('Pt', cell=[d, d, 10.], pbc=True)
	io.write(poscar, bulk, 'vasp')
>>>>>>> parent of 1d84223 (Update)

	print(f'\n{poscar} ({bulk.cell.get_bravais_lattice()}):')
	with open(poscar, 'r') as f: print(f.read())
	print('* Please add element symbol after the atomic position\n')

def gen_hspts(path, nkpts):
	bulk = io.read(poscar)

	kpts, hspts, hspts_label = bandpath(path, bulk.cell, npoints=nkpts).get_linear_kpoint_axis()
	hspts = [np.where(np.isclose(kpts, p, atol=1e-6))[0][0] for p in hspts]
	hspts_label = [l if l != 'G' else r'$\Gamma$' for l in hspts_label]

	return bulk, hspts, hspts_label

<<<<<<< HEAD
<<<<<<< HEAD
def show_hist():
	fn = f'{matdir}/{args.mat}_cell_o_HIST.nc'
	with abilab.abiopen(abidata.ref_file(fn)) as f:
		print(f'{fn}:')
		for i, (struct, etot) in enumerate(zip(f.structures, f.etotals)):
			print(f'itr: {i}  Volume: {struct.volume:.6f}  Etotal: {etot:.6f}')
		print('\n', f.final_structure, sep='', end='\n\n')

		cell = np.round(f.final_structure.lattice.abc, decimals=6)
		gen_poscar(fn='POSCAR',  cell=cell, size=(1, 1, 1))
		gen_poscar(fn='POSCAR2', cell=cell, size=(2, 2, 1))
	
=======
>>>>>>> parent of 1d84223 (Update)
=======
>>>>>>> parent of 1d84223 (Update)
def show_kpts(path, nkpts):
	bulk, hspts, hspts_label = gen_hspts(path, nkpts)

	ndivk = [hspts[i] - hspts[i-1] for i in range(1, len(hspts))]
	print('\nndivk:\n', *ndivk)

	kptbounds = get_special_points(bulk.cell)
	print('\nkptbounds:')
	for kpt in path:
		print(*kptbounds[kpt], f'#{kpt}')
	print()

def show_fatband(qnum, path, nkpts, nband):
	if len(qnum) < 2:
		l = qnum[0]
		fn = f'{matdir}/Pt_band_o_DS2_FATBANDS_at0001_Pt_is1_l000{l}'
	else:
		l, m = qnum[0], qnum[1]
		fn = f'{matdir}/Pt_band_o_DS2_FATBANDS_at0001_Pt_is1_l{l}_m{m:+}'
	print(f'Show {fn}')

	_, hspts, hspts_label = gen_hspts(path, nkpts)

	band_list = []
	with open(fn, 'r') as f:
		for line in f:
			if line.strip() and not re.match(r'[#@&]', line):
				band_list.append([float(v) for v in line.strip().split()])
	band_list = np.split(np.array(band_list), nband)

	fig, ax = plt.subplots(figsize=(8, 6), dpi=120, tight_layout=True)
<<<<<<< HEAD
<<<<<<< HEAD
	ecolor = list(mcolors.TABLEAU_COLORS)
	ylim = (-8, 5)
	for i, band_l in enumerate(band):
		for j, b in enumerate(band_l):
			ax.errorbar(b[:, 0], b[:, 1], b[:, 2], color='k', alpha=0.7, ecolor=ecolor[i], label=None if j else f'l = {i}')
			if b[0, 1] > ylim[0] and b[0, 1] < ylim[1]: ax.text(b[0, 0]+1, b[0, 1]+0.1, j+1, fontsize='xx-small')
=======
	ecolor = list(mcolors.TABLEAU_COLORS)[l if len(qnum) < 2 else m]
	for band in band_list:
		ax.errorbar(band[:, 0], band[:, 1], band[:, 2], color='k', ecolor=ecolor)
>>>>>>> parent of 1d84223 (Update)
=======
	ecolor = list(mcolors.TABLEAU_COLORS)[l if len(qnum) < 2 else m]
	for band in band_list:
		ax.errorbar(band[:, 0], band[:, 1], band[:, 2], color='k', ecolor=ecolor)
>>>>>>> parent of 1d84223 (Update)
	for p in hspts[1:-1]:
		plt.axvline(p, lw=2, color='k', alpha=0.2)
	plt.axhline(0, lw=2, ls='--', color='k', alpha=0.2)
	plt.xticks(hspts, labels=hspts_label)
	plt.xlim(0, nkpts-1)
	plt.ylim(-5, 5)
	plt.title(f'$l = {l}$' if len(qnum) < 2 else f'$l = {l}$, $m = {m}$', loc='left')
	plt.ylabel(r'$E - E_{F}$')
<<<<<<< HEAD

	fig.savefig(f'{matdir}/fig/{fn}.svg')
	plt.show()

def show_wann(path, nkpts, nband_dft, e_fermi):
	fn_wan = f'w90_band.dat'
	fn_dft = f'{args.mat}_band_o_DS2_FATBANDS_at0001_Pt_is1_l0000'

	_, hspts, hspts_label = gen_hspts(path, nkpts)

	band_wan = []
	with open(f'{matdir}/{fn_wan}', 'r') as f:
		for line in f:
			if line.strip():
				band_wan.append([float(v) for v in line.strip().split()[:2]])
	band_wan = np.array(band_wan)
	band_wan[:, 1] -= e_fermi

	x_uniq = np.unique(band_wan[:, 0])
	x_dict = {v: i for i, v in enumerate(x_uniq)}
	band_wan[:, 0] = np.array([x_dict[v] for v in band_wan[:, 0]])

	indices = np.linspace(0, len(x_uniq), nkpts).astype(int)
	res = np.array([i for i, v in enumerate(band_wan[:, 0]) if v in indices])
	band_wan = band_wan[res]

	i_dict = {v: i for i, v in enumerate(indices)}
	band_wan[:, 0] = np.array([i_dict[v] for v in band_wan[:, 0]])
	
<<<<<<< HEAD
	band_dft = []
	with open(f'{matdir}/{fn_dft}', 'r') as f:
		for line in f:
			if line.strip() and not re.match(r'[#@&]', line):
				band_dft.append([float(v) for v in line.strip().split()[:2]])
	band_dft = np.split(np.array(band_dft), nband_dft)

	fig, ax = plt.subplots(figsize=(8, 6), dpi=120, tight_layout=True)
	ax.scatter(band_wan[:, 0], band_wan[:, 1], color='r', label='wan')
	for band in band_dft: ax.plot(band[:, 0], band[:, 1], color='k')
	plt.ylim(-8, 5)
	plt.legend(fontsize='x-small', loc='upper right')

	fig.savefig(f'{matdir}/fig/wan_band.svg')
	plt.show()
		
if args.clean: clean(args.clean)
elif args.poscar: gen_poscar(fn='POSCAR', cell=[2.79, 2.79, 10.], size=(1, 1, 1))
elif args.hist: show_hist()
elif args.kpts: show_kpts(path='GXMG', nkpts=300)
elif args.wann: show_wann(path='GXMG', nkpts=300, nband_dft=40, e_fermi=-2.72)
elif args.band < 99: show_fatband(args.band, path='GXMG', nkpts=300, nband=40)
=======
=======

	fig.savefig(f'{matdir}/fig/{fn}.svg')
	plt.show()
	
>>>>>>> parent of 1d84223 (Update)
if args.clean: clean()
elif args.poscar: gen_poscar(d=2.79) #2.527276
elif args.kpts: show_kpts(path='GXMG', nkpts=100)
elif len(args.band): show_fatband(args.band, path='GXMG', nkpts=100, nband=40)
<<<<<<< HEAD
>>>>>>> parent of 1d84223 (Update)
=======
>>>>>>> parent of 1d84223 (Update)
else:
	parser.print_help()
	sys.exit(1)
