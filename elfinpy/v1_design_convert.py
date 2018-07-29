#!/usr/bin/env python3

import argparse, sys
import copy
from elfin_graph import ElfinGraph
from elfin_node import ElfinNode
from utilities import *

def compute_old_graph_txm(xdb, graph):
    nodes = graph.nodes
    doubles_data = xdb['doubles_data']
    for i in range(len(nodes)-1):
        node_a = nodes[i] 
        node_b = nodes[i+1]
        rel = doubles_data[node_a.name][node_b.name]
        for j in range(i+1):
            nodes[j].transform(rel['rot'], rel['tran'])


def parse_args(args):
    parser = argparse.ArgumentParser(description='Converts old Elfin core intermediate output into new format');
    parser.add_argument('input') # No dash means mandatory
    parser.add_argument('--output')
    parser.add_argument('--xdb_path', default='resources/xDB.json')
    parser.add_argument('--multichain_test', action='store_true')
    return parser.parse_args(args)

def main(test_args=None):
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    # Elfin core output
    ec_out = read_json(args.input)

    # Make sure we're working with the old format
    keys = ec_out.keys()
    if not 'nodes' in keys:
        print('Input file does not look like the old Elfin core output file')
        return 1

    n_nodes = len(ec_out['nodes'])
    nodes = [ 
                ElfinNode(
                i, 
                el, 
                trim=[(False if i == 0 else True), (False if i == n_nodes - 1 else True)],
                cterm_node_id=((i+1) if i < n_nodes - 1 else -1)
                ) for (i, el) in enumerate(ec_out['nodes'])
            ]

    graph = ElfinGraph('c1', nodes) # c1 for chain number 1
    graphs = [graph]

    assert(len(graphs) == 1)

    xdb = read_json(args.xdb_path)
    map(lambda i, el: compute_old_graph_txm(xdb, el), enumerate(graphs))

    if args.multichain_test:
        graphs.append(copy.deepcopy(graph))

        # Note: flipping z direction can cause problems in PyMol 
        #   visualisation (can't view as cartoon)
        graphs[0].transform([[-1,0,0],[0,-1,0],[0,0,1]],[100,100,0])
        graphs[1].transform([[1,0,0],[0,1,0],[0,0,1]],[-100,-100,0])
        graphs[1].name = 'c2'

    output_file = args.output
    if output_file == None:
        output_file = args.input.replace('.json', '.new.json')

    with open(output_file, 'wb') as ofp:
        json.dump(graphs, ofp, default=lambda o: o.__dict__)
        print('Saved to: ' + output_file)

if __name__ == '__main__':
    main()