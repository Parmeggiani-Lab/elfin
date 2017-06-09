#!/usr/bin/env python

import argparse, sys
from utils import *

def main():
	ap = argparse.ArgumentParser(description='Parse chain break residue ID from filename and output loop file');
	ap.add_argument('--input')
	ap.add_argument('--range', type=int, default=2)
	ap.add_argument('--token', default='_mc_b')

	# Exit if no argument given
	if len(sys.argv) == 1:
		ap.print_help()
		sys.exit(1)

	args = ap.parse_args()

	if args.input == '':
		print 'Must provide input'
		ap.print_help()
		sys.exit(1)

	chainBreakRId = args.input[args.input.index(args.token)+len(args.token):].replace('.pdb', '')
	print 'LOOP {} {} 0 0 1'.format(
		int(chainBreakRId) - int(args.range),
		int(chainBreakRId) + 1 + int(args.range)
	)


if __name__ == '__main__':
	main()