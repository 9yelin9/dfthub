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
parser.add_argument('--clean', action='store_true')
parser.add_argument('--cell', action='store_true')
parser.add_argument('--kpts', action='store_true')
parser.add_argument('--band', type=int, nargs='+', default=-1)
args = parser.parse_args()

def clean(keep=['Pt.sh', 'Pt_cell_.abi', 'Pt_band_.abi', 'Pt_crpa_.abi', 'fig', 'log', 'run.py']):
	for fn in set(os.listdir(os.getcwd())) - set(keep):
		fn = f'{os.getcwd()}/{fn}'
		if os.path.isfile(fn): os.remove(fn)
		else: shutil.rmtree(fn)

def init_cell(d):
	bulk = Atoms('Pt', cell=[d, d, 10.], pbc=True)
	io.write('POSCAR', bulk, 'vasp')

	print(f'POSCAR ({bulk.cell.get_bravais_lattice()}):')
	with open('POSCAR', 'r') as f: print(f.read())

def show_kpts(path, nkpts):
	bulk = 0

	kpts, hspts, hspts_label = bandpath(path, bulk.cell, npoints=nkpts).get_linear_kpoint_axis()
	hspts = [np.where(np.isclose(kpts, p, atol=1e-6))[0][0] for p in hspts]
	hspts_label = [l if l != 'G' else r'$\Gamma$' for l in hspts_label]

	ndivk = [hspts[i] - hspts[i-1] for i in range(1, len(hspts))]
	print('\nndivk:\n', *ndivk)

	kptbounds = get_special_points(bulk.cell)
	print('\nkptbounds:')
	for kpt in path:
		print(*kptbounds[kpt], f'#{kpt}')

def show_fatband(qnum, nkpts, nband=30):
	if len(qnum) < 2:
		l = qnum[0]
		fn = f'Pt_band_o_DS2_FATBANDS_at0001_Pt_is1_l000{l}'
	else:
		l, m = qnum[0], qnum[1]
		fn = f'Pt_band_o_DS2_FATBANDS_at0001_Pt_is1_l{l}_m{m:+}'
	print(f'Show {fn}')

	hspts, hspts_label = 0, 0

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

	fig.savefig(f'fig/{fn}.svg')
	plt.show()
	
if args.clean: clean()
elif args.cell: init_cell(d=2.79)
elif args.kpts: show_kpts(path='', nkpts=300)
elif len(args.band): show_fatband(args.band, nkpts=300)
else:
	parser.print_help()
	sys.exit(1)
