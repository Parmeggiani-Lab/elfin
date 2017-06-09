#!/usr/bin/env python

import glob
import Greedy
from utils import *
import time

def main():
    xDB     = readJSON('res/xDB.json')
    bmDir 	= 'bm/l10'

    designers = []
    designers.append(Greedy.GreedyDesigner(xDB, 'maxHeavy'))
    # MCDesigner
    # ...

    # Process all benchmarks
    for jsonFile in glob.glob(bmDir + '/*.json'):
    	spec = readJSON(jsonFile)
        targetLen = len(spec['nodes'])

    	for designer in designers:
            designerName = designer.__class__.__name__
            inputName = jsonFile.replace('.json', '')
            print 'Benchmarking {} on {}, target length={}'.format(
                designerName, inputName, targetLen)

            startTime = time.clock()
            (nodes,shape,score,fRot) = designer.design(inputName, spec, targetLen)
            print "{:.2f}s, score: {}".format(time.clock() - startTime, score)

            if nodes == spec['nodes']:
                print 'Pass'
            else:
                print 'Failed'
                pauseCode()
    
            makePdbFromNodes(xDB, 
                nodes, 
                'res/centered_pdb/pair', 
                jsonFile.replace('.json', suffixPdb(
                    designerName, 
                    'PCT',
                    1.0, 
                    targetLen)),
                fRot)
            
if __name__ =='__main__': safeExec(main)