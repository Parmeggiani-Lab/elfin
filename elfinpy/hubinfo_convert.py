#!/usr/bin/env python3

import argparse, sys
from collections import OrderedDict
from utilities import *


def parse_args(args):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Converts hub info metadata from csv to json.')
    parser.add_argument('input') # Absence of dash denotes mandatory argument
    return parser.parse_args(args)

def main(test_args=None):
    """main"""
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    csv_data = read_csv(args.input)
    new_data = OrderedDict({})
    for i in range(int(len(csv_data)/2)):
        row = csv_data[int(2*i)]
        chains = csv_data[int(2*i+1)]
        component_data = OrderedDict({})
        for i in range(int(len(chains)/2)):
            component_data[chains[int(2*i)]] = \
                { 'single_name': chains[int(2*i+1)], 
                'c_free': row[4]=='C_free', 
                'n_free': row[3]=='N_free' }
        new_data[row[0].replace('.pdb', '')] = \
            OrderedDict({ 'oligomer_type': row[1], 
                'symmetric': row[2] == 'symmetric', 
                'component_data': component_data})
    json.dump(new_data, 
        open(args.input.replace('.csv', '') + '.json', 'w'), 
        separators=(',',':'), 
        ensure_ascii=False, 
        indent=4)

if __name__ == '__main__':
    safe_exec(main)