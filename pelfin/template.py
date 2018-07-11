#!/usr/bin/env python3

import argparse, sys
from utilities import *

def parse_args(args):
  parser = argparse.ArgumentParser(description='Template Elfin Python script')
  parser.add_argument('input') # Absence of dash denotes mandatory argument
  return parser.parse_args(args)

def main(test_args=None):
  args = parse_args(sys.argv[1:] if test_args is None else test_args)

# def main():
#   raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
  main()