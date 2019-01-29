#!/usr/bin/env python3

import argparse, sys
import copy

try:
    from elfin_graph import ElfinGraph
    from elfin_node import ElfinNode
    from utilities import *
except ImportError as e:
    from .elfin_graph import ElfinGraph
    from .elfin_node import ElfinNode
    from .utilities import *

def compute_old_graph_txm(xdb, graph):
    nodes = graph.nodes
    double_data = xdb['double_data']
    for i in range(len(nodes)-1):
        node_a = nodes[i] 
        node_b = nodes[i+1]
        rel = double_data[node_a.name][node_b.name]
        for j in range(i+1):
            nodes[j].transform(rel['rot'], rel['tran'])


def parse_args(args):
    parser = argparse.ArgumentParser(description='Converts old Elfin core intermediate output into new format');
    parser.add_argument('input') # No dash means mandatory
    parser.add_argument('--output')
    parser.add_argument('--xdb_path', default='resources/xdb.json')
    parser.add_argument('--multichain_test', action='store_true')
    return parser.parse_args(args)

def v1_to_v2(input_json, xdb_path, multichain_test=False):# Elfin core output
    # Make sure we're working with the old format
    keys = input_json.keys()
    if not 'nodes' in keys:
        print('Error: input file is not a v1 elfin solution file.')
        exit(1)

    n_nodes = len(input_json['nodes'])
    nodes = [ 
                ElfinNode(
                    id=i, 
                    name=el, 
                    trim=[(False if i == 0 else True), (False if i == n_nodes - 1 else True)],
                    cterm_node_id=((i+1) if i < n_nodes - 1 else -1)
                ) for (i, el) in enumerate(input_json['nodes'])
            ]

    graph = ElfinGraph('c1', nodes) # c1 for chain number 1
    graphs = [graph]

    assert(len(graphs) == 1)

    xdb = read_json(xdb_path)

    for g in graphs:
        compute_old_graph_txm(xdb, g)

    if multichain_test:
        graphs.append(copy.deepcopy(graph))

        # Note: flipping z direction can cause problems in PyMol 
        #   visualisation (can't view as cartoon)
        graphs[0].transform([[-1,0,0],[0,-1,0],[0,0,1]],[100,100,0])
        graphs[1].transform([[1,0,0],[0,1,0],[0,0,1]],[-100,-100,0])
        graphs[1].name = 'c2'

    return graphs

def main(test_args=None):
    print('Deprecated. For code reference only.')
    exit()
    
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    graphs = v1_to_v2(
        read_json(args.input),
        args.xdb_path,
        args.multichain_test)

    output_file = args.output
    if output_file == None:
        output_file = args.input.replace('.json', '.v2.json')

    with open(output_file, 'w') as ofp:
        json.dump(graphs, ofp, default=lambda o: o.__dict__)
        print('Saved to: ' + output_file)

if __name__ == '__main__':
    main()