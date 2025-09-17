#!/usr/bin/env python3

import os
import sys
import time
import argparse
import numpy as np

from ase import Atoms
from ase.io import read
from ase.dft import kpoints
from ase.lattice import FCC
from ase.calculators.openmx import OpenMX

parser = argparse.ArgumentParser()
parser.add_argument('--restart', action='store_true')
args = parser.parse_args()

ngrids = 10
labels = ['L', 'G', 'X', 'W', 'L', 'K', 'G']

pts = FCC(3).get_special_points()
coords = [pts[label] for label in labels]
kpath = [[ngrids, *coords[i], *coords[i+1], labels[i], labels[i+1]] for i in range(len(coords)-1)]

bulk = read(f'{os.getcwd()}/NiO.poscar')

bulk.calc = OpenMX(
	label=f'{os.getcwd()}/openmx',
	data_path='/home/yerin/openmx3.9/DFT_DATA19',
	command='openmx',
	mpi={'processes': 8},
	definition_of_atomic_species=[['Ni', 'Ni6.0S-s2p2d1', 'Ni_PBE19S'], ['O', 'O6.0-s2p2', 'O_PBE19']],
	xc='PBE',
	maxiter=100,
	energy_cutoff=150.,
	smearing=300,
	scf_kgrid=(4, 4, 4),
	band_nkpath=len(kpath),
	band_kpath=kpath,
	convergence=1e-6,
	spinpol=True,
	#eigensolver='cluster',
	eigensolver='band',
	mixer='rmm-diis',
	#level_of_fileout=2,
	mo_fileout=True,
	scf_restart=args.restart,
)

bulk.get_potential_energy()
