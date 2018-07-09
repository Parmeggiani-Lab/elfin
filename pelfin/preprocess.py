#!/usr/bin/env python3

import argparse, sys
import subprocess, glob
from utilities import *
from pdb_utilities import *

def merge_chains(pdb):
  """
  Merge all chains in a PDB structure and re-number the residue IDs
  accordingly

  Args:
  pdb - Bio.PDB.Structure.Structure

  Returns: 
  Bio.PDB.Structure.Structure - the modified PDB, used for chaining
  calls
  """
  newChain = Bio.PDB.Chain.Chain('A')
  rid = 1
  for r in strip_residues(pdb):
    r.id = (r.id[0], rid, r.id[2])
    newChain.add(r)
    rid += 1

  # Remove old chains from model
  for model in pdb:
    ocIds = []
    for chain in model:
      ocIds.append(chain.id)
    for ocId in ocIds:
      model.detach_child(ocId)

  model.add(newChain)
  return pdb # for chaining calls

def cleanse_atoms(pdb):
  """
  Delete 1H, 2H, 3H and OXT atoms from the PDB structure

  Args:
  pdb - Bio.PDB.Structure.Structure

  Returns: 
  Bio.PDB.Structure.Structure - the modified PDB, used for chaining
  calls
  """
  for c in get_chains(pdb):
    for r in c.child_list:
      badAtoms = []
      for a in r:
        if a.name == '1H' or a.name == '2H' or a.name == '3H' or a.name == 'OXT':
          badAtoms.append(a.name)
      for ba in badAtoms:
        r.detach_child(ba)
  return pdb

def preprocess_double(double_file):
  """
  Merge chains, cleanse atoms if double_file is a simple double (no junction).
  If double_file is a complex double (one with a junction), replace its
  interfacing residues with those of a simple double, then merge chains and
  cleanse atoms.

  Loop-replacement is necessary because we have experimental data that confirm
  interface atom positions in simple doubles but not in compound doubles. If
  this is not done, the relaxation of database PDBs will result in sever
  bending and dislocation in compound doubles.

  Args:
  double_file - string path of the input double PDB file

  Returns:
  Bio.PDB.Structure.Structure - preprocessed double PDB
  """
  double_name = os.path.basename(double_file).replace('.pdb', '')
  underscores = [double_name.find('_'), double_name.rfind('_')]

  if underscores[0] == -1:
    print('Input {} is a simple double and does not need loop replacement'.format(double_file))
    return cleanse_atoms(merge_chains(read_pdb(double_file)))

  dash_idx = double_name.rfind('-')
  sdouble_first = dash_idx < underscores[0] and dash_idx < underscores[1]

  double_name_halves = [double_name[:dash_idx], double_name[dash_idx+1:]]
  sdouble_name = ''
  sdouble_name += double_name_halves[0][double_name_halves[0].rfind('_')+1:] \
        if double_name_halves[0].rfind('_') != -1 else double_name_halves[0]
  sdouble_name += '-'
  sdouble_name += double_name_halves[1][:double_name_halves[1].find('_')] \
        if double_name_halves[1].find('_') != -1 else double_name_halves[1]
  sdouble_file = double_file[:double_file.rfind('/')+1] + sdouble_name + '.pdb'

  # Load PDBs
  double = read_pdb(double_file)
  double_chains = double.child_list[0].child_list
  assert(len(double_chains) == 2)

  sdouble = read_pdb(sdouble_file) #sdouble is the simple double
  sdouble_chains = sdouble.child_list[0].child_list
  assert(len(sdouble_chains) == 2)

  # Get residue counts
  sdouble_r_count = get_pdb_residue_count(sdouble)

  # Compute interface residue range
  sdouble_start_idx = int(np.ceil(sdouble_r_count*0.375))+1 # 0.375 is 75% of first single
  sdouble_end_idx = int(np.floor(sdouble_r_count*0.625))-1 # 0.625 is 25% of second single
  sdouble_end_offset = sdouble_end_idx - len(sdouble_chains[0].child_list)

  # Find simple double middle residues
  sdouble_chain_lens = [len(c.child_list) for c in sdouble_chains]
  sdouble_mid_res = sdouble_chains[0].child_list[sdouble_start_idx:] + \
    sdouble_chains[1].child_list[:sdouble_end_offset]
  sdouble_atoms = [a for r in sdouble_mid_res for a in r.child_list if a.name == 'CA']

  # Find first half of the residues
  double_chain_lens = [len(c.child_list) for c in double_chains]
  double_start_offset = double_chain_lens[0]-(sdouble_chain_lens[0]-sdouble_start_idx)
  double_mid_res = \
    (double_chains[0].child_list[sdouble_start_idx:] if sdouble_first else \
    double_chains[0].child_list[double_start_offset:]) + \
    double_chains[1].child_list[:sdouble_end_offset]
  double_atoms = [a for r in double_mid_res for a in r.child_list if a.name == 'CA']

  # Superimpose double onto sdouble
  si = Bio.PDB.Superimposer()
  si.set_atoms(sdouble_atoms, double_atoms)
  double.transform(*si.rotran)

  # Merge chains and remove bad atoms
  cleanse_atoms(merge_chains(double))
  cleanse_atoms(merge_chains(sdouble))

  # Replace double residues where it should be sdouble residues
  double_residues = get_chain(double).child_list
  sdouble_residues = strip_residues(sdouble)
  for rIdx in range(sdouble_start_idx, sdouble_end_idx):
    offset_r_idx = rIdx + (0 if sdouble_first else double_chain_lens[0]-sdouble_chain_lens[0])
    old_r_id = double_residues[offset_r_idx].id
    double_residues[offset_r_idx] = sdouble_residues[rIdx]
    double_residues[offset_r_idx].id = old_r_id

  return double

def parse_args(args):
  parser = argparse.ArgumentParser(description='Preprocess all raw single and double PDBs');
  parser.add_argument('--input_dir', default='./resources/pdb_raw/')
  parser.add_argument('--output_dir', default='./resources/pdb_prepped/')
  parser.add_argument('--dry_run', action='store_true')
  return parser.parse_args(args)

def main(test_args=None):
  args = parse_args(sys.argv[1:] if test_args is None else test_args)

  if args.input_dir is None or args.output_dir is None:
    ap.print_help()
    sys.exit(1)

  if not args.dry_run:
    make_dir(args.output_dir)
    double_output_dir = args.output_dir + '/doubles/'
    make_dir(double_output_dir)
    single_output_dir = args.output_dir + '/singles/'
    make_dir(single_output_dir)
    hub_output_dir = args.output_dir + '/hubs/'
    make_dir(hub_output_dir)

  # Doubles
  doubleFiles = glob.glob(args.input_dir + '/doubles/*.pdb')
  N = len(doubleFiles)
  for i in range(N):
    double_file = doubleFiles[i]
    print('Prepping double [{}/{}] {}'.format(i+1, N, double_file))
    double = preprocess_double(double_file)
    if not args.dry_run:
      save_pdb(struct=double, save_path=double_output_dir + '/' + os.path.basename(double_file))

  # Singles
  singleFiles = glob.glob(args.input_dir + '/singles/*.pdb')
  N = len(singleFiles)
  for i in range(N):
    single_file = singleFiles[i]
    print('Prepping single [{}/{}] {}'.format(i+1, N, single_file))
    # Singles need nothing other than cleansing
    single = cleanse_atoms(read_pdb(single_file))
    if not args.dry_run:
      save_pdb(struct=single, save_path=single_output_dir + '/' + os.path.basename(single_file))

  # Hubs
  hubFiles = glob.glob(args.input_dir + '/hubs/*.pdb')
  N = len(hubFiles)
  for i in range(N):
    hub_file = hubFiles[i]
    print('Prepping hub [{}/{}] {}'.format(i+1, N, hub_file))
    hub = cleanse_atoms(read_pdb(hub_file))
    if not args.dry_run:
      save_pdb(struct=hub, save_path=hub_output_dir + '/' + os.path.basename(hub_file))

if __name__ == '__main__':
  safe_exec(main)