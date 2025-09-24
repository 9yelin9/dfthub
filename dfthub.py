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
parser.add_argument('--clean', action='store_true')
parser.add_argument('--poscar', action='store_true')
parser.add_argument('--kpts', action='store_true')
parser.add_argument('--band', type=int, nargs='+', default=-1)
args = parser.parse_args()

matdir = f'{os.getcwd()}/{args.mat}'
poscar = f'{matdir}/POSCAR'

def clean(keep=r'POSCAR|.*\.abi|fig'):
	for fn in os.listdir(matdir):
		fn = f'{matdir}/{fn}'
		if not re.search(keep, fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

def gen_poscar(d):
	bulk = Atoms('Pt', cell=[d, d, 10.], pbc=True)
	io.write(poscar, bulk, 'vasp')

	print(f'\n{poscar} ({bulk.cell.get_bravais_lattice()}):')
	with open(poscar, 'r') as f: print(f.read())
	print('* Please add element symbol after the atomic position\n')

def gen_hspts(path, nkpts):
	bulk = io.read(poscar)

	kpts, hspts, hspts_label = bandpath(path, bulk.cell, npoints=nkpts).get_linear_kpoint_axis()
	hspts = [np.where(np.isclose(kpts, p, atol=1e-6))[0][0] for p in hspts]
	hspts_label = [l if l != 'G' else r'$\Gamma$' for l in hspts_label]

	return bulk, hspts, hspts_label

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
	ecolor = list(mcolors.TABLEAU_COLORS)[l if len(qnum) < 2 else m]
	for band in band_list:
		ax.errorbar(band[:, 0], band[:, 1], band[:, 2], color='k', ecolor=ecolor)
	for p in hspts[1:-1]:
		plt.axvline(p, lw=2, color='k', alpha=0.2)
	plt.axhline(0, lw=2, ls='--', color='k', alpha=0.2)
	plt.xticks(hspts, labels=hspts_label)
	plt.xlim(0, nkpts-1)
	plt.ylim(-5, 5)
	plt.title(f'$l = {l}$' if len(qnum) < 2 else f'$l = {l}$, $m = {m}$', loc='left')
	plt.ylabel(r'$E - E_{F}$')

	fig.savefig(f'{matdir}/fig/{fn}.svg')
	plt.show()
	
if args.clean: clean()
elif args.poscar: gen_poscar(d=2.79) #2.527276
elif args.kpts: show_kpts(path='GXMG', nkpts=100)
elif len(args.band): show_fatband(args.band, path='GXMG', nkpts=100, nband=40)
else:
	parser.print_help()
	sys.exit(1)
