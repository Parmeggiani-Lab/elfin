#!/usr/bin/env python

import argparse, sys
import numpy as np
import Kabsch
from utils import *

def getSpecSolRot(specFile, solCSV):
    if specFile.rfind('.csv') != -1:
        specPts = readCSVPoints(specFile)
    elif specFile.rfind('.json') != -1:
        with open(specFile, 'r') as file:
            specPts = np.asarray(json.load(file)['coms'])
    else:
        print 'Unknown spec file format'

    solPts = readCSVPoints(solCSV)

    # Centre both pts to their ends
    centredSpec = specPts - specPts[-1]
    centredSol = solPts - solPts[-1]

    # Equalise sample points
    solUpPts = upsample(centredSpec, centredSol)
    solUpPts = solUpPts - solUpPts[-1]

    # Find Kabsch rotation for solution -> spec
    rot = Kabsch.kabsch(solUpPts, centredSpec)

    rotTp = np.transpose(rot)
    rotTpTran = np.append(rotTp, np.transpose([[0,0,0]]), axis=1)
    pymolRotMat = np.append(rotTpTran, [[0,0,0,1]], axis=0)
    pymolRotMatStr = '[' + ', '.join(map(str, pymolRotMat.ravel())) + ']'

    return pymolRotMatStr

def main():
	ap = argparse.ArgumentParser(description='Generate Grid Search configurations');
	ap.add_argument('specFile')
	ap.add_argument('solFile')
	args = ap.parse_args()

	print getSpecSolRot(args.specFile, args.solFile)

if __name__ == '__main__':
	main()