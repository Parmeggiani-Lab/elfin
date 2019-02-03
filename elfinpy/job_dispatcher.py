#!/usr/bin/env python3

import argparse, sys
import multiprocessing
import subprocess

from utilities import *

def dispatch(*cmd_and_arg):
    """Dispatches a process to run cmd with given arguments."""
    subprocess.check_call(*cmd_and_arg)

def parse_args(args):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Run single threaded jobs in separate processes')
    parser.add_argument('cmd_list')
    parser.add_argument('-worker_count', default='cpu_count')
    return parser.parse_args(args)

def main(test_args=None):
    """main"""
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    try:
        if args.worker_count == 'cpu_count':
            pool = multiprocessing.Pool(multiprocessing.cpu_count())
        else:
            pool = multiprocessing.Pool(int(args.worker_count))
    except Exception as e:
        raise e

    pool.map(dispatch, read_csv(args.cmd_list, delim=' '))

if __name__ == '__main__':
    safe_exec(main)