#!/usr/bin/env python3

import argparse, sys
import multiprocessing
import subprocess

from utilities import *

def star_call(*cmd_and_arg):
  subprocess.check_call(*cmd_and_arg)

def main():
  ap = argparse.ArgumentParser(description='Run single threaded jobs in separate processes')
  ap.add_argument('cmd_list')
  ap.add_argument('-worker_count', default='cpu_count')
  args = ap.parse_args()

  try:
    if args.worker_count == 'cpu_count':
      pool = multiprocessing.Pool(multiprocessing.cpu_count())
    else:
      pool = multiprocessing.Pool(int(args.worker_count))
  except Exception as e:
    raise e

  pool.map(star_call, read_csv(args.cmd_list, delim=' '))

if __name__ == '__main__':
  safe_exec(main)