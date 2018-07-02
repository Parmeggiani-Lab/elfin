#!/usr/bin/env python3

import argparse, sys
from utilities import *

def main():
  ap = argparse.ArgumentParser(description='Template Elfin Python script')
  ap.add_argument('input') # Absence of dash denotes mandatory argument
  args = ap.parse_args()

if __name__ == '__main__':
  safe_exec(main)