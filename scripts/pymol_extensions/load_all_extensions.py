#!/usr/bin/env python2

#
# This is a PyMol extension script to load all PyMol extensions in the same
# directory
#

def main():
  raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
  main()

def __in_pymol():
  import sys, os

  curr_dir = os.getcwd()
  sys.path.append(os.path.join(curr_dir, os.pardir, os.pardir)) # for pelfin

  import importlib
  ext_names = [
    'draw_lines',
    'batch_convert',
    'transform_helper'
  ]

  for ext in ext_names:
    if ext in sys.modules:
      reload(sys.modules[ext])
    else:
      importlib.import_module(ext)

  print('All Extensions Loaded')

in_pymol = False
try:
  import pymol
  in_pymol = True
except ImportError as ie:
  main()

if in_pymol: __in_pymol()