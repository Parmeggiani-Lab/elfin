#!/usr/bin/env python

import argparse, sys
import numpy as np
import kabsch
from utilities import *

def get_spec_sol_rot(spec_file, sol_csv):
    if spec_file.rfind('.csv') != -1:
        spec_pts = read_csv_points(spec_file)
    elif spec_file.rfind('.json') != -1:
        with open(spec_file, 'r') as file:
            spec_pts = np.asarray(json.load(file)['coms'])
    else:
        print 'Unknown spec file format'

    solPts = read_csv_points(sol_csv)

    # Centre both pts to their ends
    centred_spec = spec_pts - spec_pts[-1]
    centred_sol = solPts - solPts[-1]

    # Equalise sample points
    sol_up_pts = upsample(centred_spec, centred_sol)
    sol_up_pts = sol_up_pts - sol_up_pts[-1]

    # Find Kabsch rotation for solution -> spec
    rot = Kabsch.kabsch(sol_up_pts, centred_spec)
    return gen_pymol_txm(rot)

def main():
	ap = argparse.ArgumentParser(
        description='Generate spec to solution rotation string for Pymol')
	ap.add_argument('spec_file')
	ap.add_argument('sol_file')
	args = ap.parse_args()

	print(get_spec_sol_rot(args.spec_file, args.sol_file))

if __name__ == '__main__':
	main()