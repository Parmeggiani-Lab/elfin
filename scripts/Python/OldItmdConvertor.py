#!/usr/bin/env python

import argparse, sys
import copy
from ElfinUtils import *

def computeOldGraphTxm(xdb, graph):
    nodes = graph.nodes
    pairsData = xdb['pairsData']
    for i in xrange(0, len(nodes)-1):
        nodeA = nodes[i] 
        nodeB = nodes[i+1]
        rel = pairsData[nodeA.name][nodeB.name]
        for j in xrange(0, i+1):
            nodes[j].transform(rel['rot'], rel['tran'])

def main():
    ap = argparse.ArgumentParser(description='Converts old Elfin core intermediate output into new format');
    ap.add_argument('input') # No dash means mandatory
    ap.add_argument('--output')
    ap.add_argument('--xdbPath', default='res/xDB.json')
    ap.add_argument('--multichainTest', action='store_true')
    args = ap.parse_args()
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)

    ecOut = readJSON(args.input)

    # Make sure we're working with the old format
    keys = ecOut.keys()
    if not 'nodes' in keys:
        print('Input file does not look like the old Elfin core output file')
        return 1

    nNodes = len(ecOut['nodes'])
    nodes = map(
        lambda (i, el): ElfinNode(
            i, 
            el, 
            trim=[(False if i == 0 else True), (False if i == nNodes - 1 else True)],
            ctermNodeId=((i+1) if i < nNodes - 1 else -1)
            ), 
        enumerate(ecOut['nodes'])
    )

    graph = ElfinGraph('c1', nodes) # c1 for chain number 1
    graphs = [graph]

    assert(len(graphs) == 1)

    xdb = readJSON(args.xdbPath)
    map(lambda (i, el): computeOldGraphTxm(xdb, el), enumerate(graphs))

    if args.multichainTest:
        graphs.append(copy.deepcopy(graph))

        # Note: flipping z direction can cause problems in PyMol 
        #   visualisation (can't view as cartoon)
        graphs[0].transform([[-1,0,0],[0,-1,0],[0,0,1]],[100,100,0])
        graphs[1].transform([[1,0,0],[0,1,0],[0,0,1]],[-100,-100,0])
        graphs[1].name = 'c2'

    outputFile = args.output
    if outputFile == None:
        outputFile = args.input.replace('.json', '.new.json')

    with open(outputFile, 'wb') as ofp:
        json.dump(graphs, ofp, default=lambda o: o.__dict__)
        print 'Saved to: ' + outputFile

if __name__ == '__main__':
    main()