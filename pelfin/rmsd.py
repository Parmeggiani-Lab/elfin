#!/usr/bin/env python3

import argparse, sys
from utilities import *

def main():
  ap = argparse.ArgumentParser(description='Sliding-window RMSD calculator')
  ap.add_argument('solution_dir') # Absence of dash denotes mandatory argument
  ap.add_argument('minimised_dir')
  ap.add_argument('-window_len', default=300)
  ap.add_argument('-overlap_ratio', default=0.5)
  ap.add_argument('-warn_threshold', default=5.0)
  args = ap.parse_args()

  if(args.window_len < 0)
    raise ValueError('Invalid window_len: must be greater than 0')

  if(args.overlap_ratio < 0 or args.overlap_ratio > 1.0)
    raise ValueError('Invalid overlap_ratio: must be between 0 and 1.0')

  if(args.warn_threshold < 0)
    raise ValueError('Invalid warn_threshold: must be greater than 0')

  overlap = int(args.overlap_ratio * args.window_len)
  print('Using window length {} and overlap ratio {} (={} CA atoms), warn at {}A'.format(
    args.window_len, args.overlap_ratio, overlap, args.warn_threshold))

  minimised_files = glob.glob(args.minimised_dir + '/*.pdb')

  for i in range(0, len(minimised_files)):
    minimised_file = minimised_files[i]
    solution_file = args.solution_dir + minimised_file[minimised_file.rfind('/'):].replace('_0001.pdb', '.pdb')
      
    minimised_pdb = read_pdb(minimised_file)
    solution_pdb = read_pdb(solution_file)

    minimised_CAs = []
    solution_CAs = []
    for a in minimised_pdb.get_atoms():
      if a.id == 'CA':
        minimised_CAs.append(a)
    for a in solution_pdb.get_atoms():
      if a.id == 'CA':
        solution_CAs.append(a)

    # Superimpose the two structures before comparing
    si = Bio.PDB.Superimposer()
    si.set_atoms(minimised_CAs, solution_CAs)
    rot, tran = si.rotran

    solution_pdb.transform(rot, tran)

    minimised_CA_coords = []
    solution_CA_coords = []
    for a in minimised_pdb.get_atoms():
      if a.id == 'CA':
        minimised_CA_coords.append(a.coord)
    for a in solution_pdb.get_atoms():
      if a.id == 'CA':
        solution_CA_coords.append(a.coord)

    n_minimised_CAs =len(minimised_CA_coords)

    if n_minimised_CAs != len(solution_CA_coords):
      print('{}: Fatal! Number of CAs are different... Solution: {}, Minimised: {}'.format(
        solution_file, len(solution_CA_coords), n_minimised_CAs))
    else:
      stats = {
        'min': float('inf'),
        'avg': float('nan'),
        'max': -float('inf')
      }

      sum_win_rmsd = 0.0
      start_index = 0
      window_count = 0
      while start_index + args.window_len - 1 < n_minimised_CAs:
        sumD = 0.0
        for i in range(start_index, start_index+args.window_len):
          sumD += np.linalg.norm(minimised_CA_coords[i] - solution_CA_coords[i], 2)

        winRmsd = np.sqrt(sumD / args.window_len)
        if winRmsd < stats['min']:
          stats['min'] = winRmsd
        if winRmsd > stats['max']:
          stats['max'] = winRmsd

        sum_win_rmsd += winRmsd
        window_count += 1

        start_index += overlap

      stats['avg'] = sum_win_rmsd / window_count

      print('{:20} avg: {:10} min {:10} max {:10}'.format(
        solution_file, stats['avg'], stats['min'], stats['max']))

      if stats['max'] > args.warn_threshold:
        print('Warning: max window RMSD of {} exceeded!'.format(args.warn_threshold))

if __name__ == '__main__':
  safe_exec(main)
