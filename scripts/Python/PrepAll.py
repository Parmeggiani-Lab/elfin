#!/usr/bin/env python

import argparse, sys
import subprocess, glob
from shutil import copyfile
from ElfinUtils import *

def main():
	ap = argparse.ArgumentParser(description='Preprocess all raw single and double PDBs');
	ap.add_argument('--inputDir', default='./resources/pdb_raw/')
	ap.add_argument('--outputDir', default='./resources/pdb_prepped/')
	args = ap.parse_args()

	if args.inputDir is None or args.outputDir is None:
		ap.print_help()
		sys.exit(1)

	mkdir(args.outputDir)
	doubleOutputDir = args.outputDir + '/doubles/'
	mkdir(doubleOutputDir)
	singleOutputDir = args.outputDir + '/singles/'
	mkdir(singleOutputDir)

	# Assume inputDir has single and double folders
	doubleFiles = glob.glob(args.inputDir + '/doubles/*.pdb')
	singleFiles = glob.glob(args.inputDir + '/singles/*.pdb')

	# Pre-process (replace interface) each double PDB
	N = len(doubleFiles)
	for i in range(N):
		print 'Prepping [{}/{}] {}'.format(i+1, N, doubleFiles[i])
		subprocess.check_call(
			[
				'./scripts/Python/PrepDouble.py', 
				'--outputDir', 
				doubleOutputDir, 
				doubleFiles[i]
			])

	# Simply copy singles over since they do not need to be modified
	for sf in singleFiles:
		copyfile(sf, singleOutputDir + sf.split('/')[-1])

if __name__ == '__main__':
	main()