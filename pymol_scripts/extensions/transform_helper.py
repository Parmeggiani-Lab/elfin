#!/usr/bin/env python3

#
# A PyMol extension script to shorten the transform_selection() command
#

def main():
  raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
  main()

in_pymol = False
try:
  import pymol
  in_pymol = True
except ImportError as ie:
  main()

if in_pymol: 
  from pymol import cmd

  import numpy as np

  @cmd.extend
  def tx(
      obj_name=None, 
      tran_before=[0,0,0], 
      rot=[[1,0,0],[0,1,0],[0,0,1]], 
      tran_after=[0,0,0]
    ):
    '''
    Transforms an object.

    Args:
    - obj_name - string
    - tran_before - a 3x1 translation vector applied before rotation
    - rot - a 3x3 rotation matrix
    - tran_after - a 3x1 translation vector applied after rotation
    '''
    if obj_name is None:
      print(tx.__doc__)
    else:
      rot_tran_mat = np.transpose(rot)
      rot_tran_mat = np.append(rot_tran_mat, np.transpose([tran_after]), axis=1)
      rot_tran_mat = np.append(rot_tran_mat, [tran_before + [1]], axis=0)
      pymol_rot_tran_vals = [v for row in rot_tran_mat for v in row]
      cmd.transform_selection(obj_name, matrix=pymol_rot_tran_vals, homogenous=0)

  @cmd.extend
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

  print('Transform Helper Loaded')