#!/usr/bin/env python

import argparse, sys
import numpy as np
from ElfinUtils import *

def main():
	ap = argparse.ArgumentParser(description='Generate Grid Search configurations');
	ap.add_argument('pairName')
	ap.add_argument('--xdbPath', default='res/xDB.json')
	ap.add_argument('--pairDir', default='/Users/joy/src/elfin/res/aligned_modules/pair/')

	args = ap.parse_args()

	xDB = readJSON(args.xdbPath)

	resetString = \
		('delete {}\n' + \
		'load {}\n' + \
		'hide everything, {}\n' + \
		'show cartoon, {}\n') \
		.format(
			args.pairName, 
			args.pairDir + '/' + args.pairName + '.pdb',
			args.pairName,
			args.pairName
		)
	singleNames = args.pairName.split('-')
	rel = xDB['pairsData'][singleNames[0]][singleNames[1]]
	rotTp = np.transpose(rel['rot'])
	rotTpTran = np.append(rotTp, np.transpose([rel['tran']]), axis=1)
	pymolRotMat = np.append(rotTpTran, [[0,0,0,1]], axis=0)

	pymolRotMatStr = '[' + ', '.join(map(str, pymolRotMat.ravel())) + ']'

	txString = \
		'cmd.transform_selection({}, {}, homogenous=0)' \
		.format(
			"\'" + args.pairName + "\'",
			pymolRotMatStr
		)
	
	print resetString
	print txString

if __name__ == '__main__':
	main()