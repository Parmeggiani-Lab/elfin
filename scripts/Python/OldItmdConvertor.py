#!/usr/bin/env python

import argparse, sys
from ElfinUtils import *

def computeOldGraphTxm(xdb, graph):
    nodes = graph.nodes
    pairsData = xdb['pairsData']
    for i in xrange(0, len(nodes) - 1):
        nodeA = nodes[i] 
        nodeB = nodes[i+1]
        rel = pairsData[nodeA.name][nodeB.name]

        for j in xrange(0, i + 1):
            nodes[j].rot = np.dot(nodes[j].rot, rel['rot'])
            nodes[j].tran = np.dot(nodes[j].tran, rel['rot']) + rel['tran']

    # Convert np array back to normal list for JSON stringifying
    for i in xrange(0, len(nodes) - 1):
        n = nodes[i]
        n.rot = n.rot.tolist()
        n.tran = n.tran.tolist()

def main():
    ap = argparse.ArgumentParser(description='Converts old Elfin core intermediate output into new format');
    ap.add_argument('input') # No dash means mandatory
    ap.add_argument('--output')
    ap.add_argument('--xdbPath', default='res/xDB.json')
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
        lambda(i,el): ElfinNode(
            i, 
            el, 
            isStart=(True if i == 0 else False),
            ctermNodes=([i+1] if i < nNodes - 1 else [])
            ), 
        enumerate(ecOut['nodes'])
    )

    graph = ElfinGraph('c1', nodes, nodes[0].id) # c1 for chain number 1
    graphs = [graph]

    assert(len(graphs) == 1)

    xdb = readJSON(args.xdbPath)
    map(lambda(i, el): computeOldGraphTxm(xdb, el), enumerate(graphs))

    for i in xrange(0, len(graphs[0].nodes)):
        n = nodes[i]
        print 'Node #' + str(i) + ' ' + genPymolTxm(n.rot, n.tran)

    # pauseCode()
    outputFile = args.output
    if outputFile == None:
        outputFile = args.input.replace('.json', '.new.json')

    with open(outputFile, 'wb') as ofp:
        json.dump(graphs, ofp, default=lambda(o): o.__dict__)
        print 'Saved to: ' + outputFile

if __name__ == '__main__':
    main()