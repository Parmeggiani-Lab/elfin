#!/usr/bin/env python

import argparse, sys
import subprocess, glob
from shutil import copyfile
from utils import *

def main():
	ap = argparse.ArgumentParser(description='Template Python script');
	ap.add_argument('--inputDir', default='./res/db/')
	ap.add_argument('--outputDir', default='./res/preprocessed/')
	args = ap.parse_args()

	if args.inputDir is None or args.outputDir is None:
		ap.print_help()
		sys.exit(1)

	mkdir(args.outputDir)
	pairOutputDir = args.outputDir + '/pair/'
	mkdir(pairOutputDir)
	singleOutputDir = args.outputDir + '/single/'
	mkdir(singleOutputDir)

	# Assume inputDir has single and pair folders
	pairFiles = glob.glob(args.inputDir + '/pair/*.pdb')
	singleFiles = glob.glob(args.inputDir + '/single/*.pdb')

	# Pre-process (replace interface) each pair PDB
	for pf in pairFiles:
		subprocess.check_call(
			[
				'./scripts/Python/PreprocessPairPDB.py', 
				'--outputDir', 
				pairOutputDir, 
				pf
			])

	# Simply copy singles over since they do not need to be modified
	for sf in singleFiles:
		copyfile(sf, singleOutputDir + sf[sf.rfind('/')+1:])

if __name__ == '__main__':
	main()