#!/usr/bin/env python

import glob
from Designer import *
import Greedy
from utils import *
import time
import numpy as np

def main():
    xDB     = readJSON('resources/xDB.json')
    bmDir 	= 'bm/l10'

    designers = []
    designers.append(Greedy.GreedyDesigner(xDB, 'maxHeavy'))
    # MCDesigner
    # ...

    # Process all benchmarks
    for jsonFile in glob.glob(bmDir + '/*.json'):
    	spec = readJSON(jsonFile)

        # we add one point between each pair of com
        coms = np.asarray(spec['coms'])
        targetLen = len(coms)
        mids = np.asarray([np.mean(p, axis=0) for p in zip(coms, np.roll(coms, -1, axis=0))[0:-1]])

        spec['coms'] = np.append([val for p in zip(coms, mids) for val in p], [coms[-1]], axis=0)

    	for designer in designers:
            designerName = designer.__class__.__name__
            print 'Benchmarking {} on {}, target length={}'.format(
                designerName, jsonFile, targetLen)

            startTime = time.clock()
            (nodes,shape,score) = designer.design(spec, targetLen)
            print "{:.2f}s, score: {}".format(time.clock() - startTime, score)

            if nodes == spec['nodes']:
                print 'Pass'
            else:
                print 'Failed'
                # pauseCode()

            makePdbFromNodes(xDB, nodes, 'resources/centered_pdb/pair', jsonFile.replace('.json', '_' + designerName + '_LV2.pdb'))

if __name__ =='__main__': safeExec(main)