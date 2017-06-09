#!/usr/bin/env python

import argparse, sys
from utils import *

def main():
	ap = argparse.ArgumentParser(description='Template Python script');
	ap.add_argument('input') # No dash means mandatory
	args = ap.parse_args()
	
	pdb = readPdb('whatever', args.input)
	print '{} has {} residues'.format(args.input, getResidueCount(pdb))

if __name__ == '__main__':
	main()