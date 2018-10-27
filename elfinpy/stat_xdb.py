#!/usr/bin/env python3

import argparse, sys
from utilities import *

def parse_args(args):
        parser = argparse.ArgumentParser(description='Prints module radii stat from xdb')
        return parser.parse_args(args)

def main(test_args=None):
        args = parse_args(sys.argv[1:] if test_args is None else test_args)
        (avg_d, min_d, max_d) = com_dist_info(read_json('resources/xdb.json'))
        print('Distances avg: {}, min: {}, max: {}'.format(avgD, minD, maxD))

if __name__ =='__main__': 
        safe_exec(main)