#!/usr/bin/env python3

#
# A PyMol extension script to shorten the transform_selection() command
#

def main():
  raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
  main()

def __in_pymol():
  from pymol import cmd

  import numpy as np

  def tx(obj_name=None, rot=[[1,0,0],[0,1,0],[0,0,1]], tran=[0,0,0]):
    '''
    Transforms an object.

    Args:
    - obj_name - string
    - rot - a 3x3 rotation matrix
    - tran - a 3x1 translation vector
    '''
    if obj_name is None:
      print(tx.__doc__)
    else:
      rot_tp = np.transpose(rot)
      rot_tp_tran = np.append(rot_tp, np.transpose([tran]), axis=1)
      pymol_rot_mat = np.append(rot_tp_tran, [[0,0,0,1]], axis=0)
      cmd.transform_selection(obj_name, matrix=pymol_rot_mat, homogenous=0)

  def multi_tx(obj_name=None, rottran_list=[]):
    '''
    Transforms an object over a list of rot, tran tuples.

    Args:
    - obj_name - string
    - rottran_list - a list of (rot, tran) lists or tuples
    '''
    if obj_name is None:
      print(multi_tx.__doc__)
    else:
      for rot, tran in rottran_list:
        etc(obj_name, rot=rot, tran=tran)

  cmd.extend("tx", tx)
  cmd.extend("multi_tx", multi_tx)

  print('Transform Helper Loaded')

in_pymol = False
try:
  import pymol
  in_pymol = True
except ImportError as ie:
  main()

if in_pymol: __in_pymol()