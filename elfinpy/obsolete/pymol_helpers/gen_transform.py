#!/usr/bin/env python

import argparse, sys
import numpy as np
from utilities import *

def main():
	ap = argparse.ArgumentParser(description='Generate Pymol Transform');
	ap.add_argument('double_name')
	ap.add_argument('--xdb_path', default='./resources/xdb.json')
	ap.add_argument('--double_dir', default='./resources/pdb_aligned/doubles/')

	args = ap.parse_args()

	xDB = readJSON(args.xdb_path)

	reset_string = \
		('delete {}\n' + \
		'load {}\n' + \
		'hide everything, {}\n' + \
		'show cartoon, {}\n') \
		.format(
			args.double_name, 
			args.double_dir + '/' + args.double_name + '.pdb',
			args.double_name,
			args.double_name
		)
	single_names = args.double_name.split('-')
	rel = xDB['doublesData'][single_names[0]][single_names[1]]
	pymol_rot_mat_str = gen_pymol_txm(rel['rot'], rel['tran'])

	tx_string = \
		'cmd.transform_selection({}, {}, homogenous=0)' \
		.format(
			"\'" + args.double_name + "\'",
			pymol_rot_mat_str
		)
	
	print(reset_string)
	print(tx_string)

if __name__ == '__main__':
	main()