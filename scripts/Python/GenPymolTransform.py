#!/usr/bin/env python

import argparse, sys
import numpy as np
from ElfinUtils import *

def main():
	ap = argparse.ArgumentParser(description='Generate Pymol Transform');
	ap.add_argument('doubleName')
	ap.add_argument('--xdbPath', default='resources/xDB.json')
	ap.add_argument('--doubleDir', default='/Users/joy/src/elfin/resources/pdb_aligned/doubles/')

	args = ap.parse_args()

	xDB = readJSON(args.xdbPath)

	resetString = \
		('delete {}\n' + \
		'load {}\n' + \
		'hide everything, {}\n' + \
		'show cartoon, {}\n') \
		.format(
			args.doubleName, 
			args.doubleDir + '/' + args.doubleName + '.pdb',
			args.doubleName,
			args.doubleName
		)
	singleNames = args.doubleName.split('-')
	rel = xDB['doublesData'][singleNames[0]][singleNames[1]]
	pymolRotMatStr = genPymolTxm(rel['rot'], rel['tran'])

	txString = \
		'cmd.transform_selection({}, {}, homogenous=0)' \
		.format(
			"\'" + args.doubleName + "\'",
			pymolRotMatStr
		)
	
	print resetString
	print txString

if __name__ == '__main__':
	main()