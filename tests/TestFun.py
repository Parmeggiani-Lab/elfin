#!/usr/bin/env python

import glob
from Designer import *
import Greedy
from utils import *
import numpy as np
import re
import time

def main():
    xDB     = readJSON('resources/xDB.json')
    bmDir 	= 'bm/fun'

    designers = []
    designers.append(Greedy.GreedyDesigner(xDB, 'maxHeavy'))
    # MCDesigner
    # ...

    (avgD, minD, maxD) = getXDBStat()

    # Process all benchmarks
    for charFile in glob.glob(bmDir + '/*.csv'):
        with open(charFile, 'r') as file:
    	    pts = [[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n')]

        spec = {'coms': pts}
        npts = np.asarray(pts)
        # Use total length/avgD as heuristic. avgD is average xDB pair distance
        targetLen = int(np.ceil(sum(np.linalg.norm(npts-np.roll(npts, 1, axis=0), axis=1)) / avgD))

    	for designer in designers:
            designerName = designer.__class__.__name__
            print 'Benchmarking {} on {}, target length={}'.format(
                designerName, charFile, targetLen)

            startTime = time.clock()
            (nodes,shape,score) = designer.design(spec, targetLen)
            print "{:.2f}s, score: {}".format(time.clock() - startTime, score)
            
            makePdbFromNodes(xDB, nodes, 'resources/centered_pdb/pair', charFile.replace('.csv', '_' + designerName + '_FUN.pdb'))
            
if __name__ =='__main__': safeExec(main)