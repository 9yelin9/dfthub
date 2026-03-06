#!/usr/bin/env python3

import os
import re
import sys
import time
import shutil
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--clean', type=str, nargs=2)
args = parser.parse_args()

def clean(target_dir, target_type):
	for fn in [f'{target_dir}/{fn}' for fn in os.listdir(target_dir) if re.search(f'{target_type}_', fn)]:
		if not re.search('[.]abi|[.]win', fn):
			if os.path.isfile(fn): os.remove(fn)
			else: shutil.rmtree(fn)

if args.clean: clean(*args.clean)
else:
	parser.print_help()
	sys.exit(1)
