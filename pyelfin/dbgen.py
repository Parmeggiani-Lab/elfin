#!/usr/bin/env python

import Bio.PDB
import glob
import numpy as np
import codecs, json
import argparse
from collections import OrderedDict

from utilities import *

def main():
  ap = argparse.ArgumentParser(description='Generates the xDB database from preprocessed single and double modules.');
  ap.add_argument('--relaxed_pdbs_dir', default='./resources/pdb_relaxed/')
  ap.add_argument('--output', default='./resources/xdb.json')
  ap.add_argument('--aligned_pdb_dir', default='./resources/pdb_aligned/')
  args = ap.parse_args()

  XDBGenerator(
    args.relaxed_pdbs_dir, 
    args.aligned_pdb_dir, 
    args.output
  ).run()

class XDBGenerator:
  def __init__(
    self,
    relaxed_pdbs_dir,
    aligned_pdb_dir,
    out_file
  ):
    self.relaxed_pdbs_dir = relaxed_pdbs_dir
    make_dir(aligned_pdb_dir)
    make_dir(aligned_pdb_dir + '/doubles/')
    make_dir(aligned_pdb_dir + '/singles/')
    make_dir(aligned_pdb_dir + '/hubs/')
    self.aligned_pdb_dir  = aligned_pdb_dir
    self.out_file          = out_file
    self.si               = Bio.PDB.Superimposer()
    self.doubles_data     = {}
    self.singles_data     = {}
    self.single_pdbs      = {}

  def get_centre_of_mass(
    self, 
    child, 
    mother=None, 
    child_resi_offset=0, 
    mother_resi_offset=0,
    match_count=-1
  ):
    CAs = []
    for a in child.get_atoms():
      if(a.name == 'CA'):
        CAs.append(a.get_coord().astype('float64'))
    com = np.mean(CAs, axis=0)

    if mother is not None:
      # This is for finding COM of a single inside a double
      _, tran = self.get_rot_trans(
        child, 
        mother, 
        moving_resi_offset=child_resi_offset, 
        fixed_resi_offset=mother_resi_offset,
        match_count=match_count
      )

      com += tran
    return com

  def move_to_origin(self, pdb):
    com = self.get_centre_of_mass(pdb)

    # No rotation - just move to centre
    pdb.transform([[1,0,0],[0,1,0],[0,0,1]], -com)

  def align(
    self, 
    moving, 
    fixed, 
    moving_resi_offset=0, 
    fixed_resi_offset=0,
    match_count=-1
  ):
    rot, tran = self.get_rot_trans(
      moving, 
      fixed,
      moving_resi_offset=moving_resi_offset,
      fixed_resi_offset=fixed_resi_offset,
      match_count=match_count
    )
    moving.transform(rot, tran)

  def get_rot_trans(
    self, 
    moving, 
    fixed, 
    moving_resi_offset=0, 
    fixed_resi_offset=0,
    match_count=-1
  ):
    moving_chain = get_chain(moving)
    ma = [
          al[0] for al in [[a for a in r.child_list if a.name == 'CA'] 
          for r in moving_chain.child_list[moving_resi_offset:(moving_resi_offset+match_count)]]
        ]

    fixed_chain = get_chain(fixed)
    fa = [
          al[0] for al in [[a for a in r.child_list if a.name == 'CA'] 
          for r in fixed_chain.child_list[fixed_resi_offset:(fixed_resi_offset+match_count)]]
        ]

    self.si.set_atoms(fa, ma)

    #   The rotation from BioPython is the second dot operand instead of
    # the conventional first dot operand.
    #
    #   This means instead of R*v + T, the actual transform is done with
    # v'*R + T
    #
    #   This is important to understand why I did the rotation maths this
    # way in the C++ GA
    return self.si.rotran

  def get_radii(self, pose):
    # Warning: this function assumes pose is centered!

    natoms = 0;
    rgSum = 0;
    max_ca_dist = 0;

    nHeavy = 0;
    max_heavy_dist = 0;
    for a in pose.get_atoms():
      dist = np.linalg.norm(
        a.get_coord().astype('float64'));

      rgSum += dist;

      if(a.name =='CA'):
        max_ca_dist = max(max_ca_dist, dist);

      if(a.element != 'H'):
        max_heavy_dist = max(max_heavy_dist, dist);
        nHeavy = nHeavy + 1;

      natoms = natoms + 1;

    average_all = rgSum / natoms;
    return OrderedDict([
      ('average_all', average_all),
      ('max_ca_dist', max_ca_dist),
      ('max_heavy_dist', max_heavy_dist)
    ]);

  def process_single(self, file_name):
    single_name = file_name.split('/')[-1].replace('.pdb', '')
    single = read_pdb(file_name)
    self.move_to_origin(single)
    save_pdb(
      single, 
      self.aligned_pdb_dir + '/singles/' + single_name + '.pdb'
    )
    self.single_pdbs[single_name] = single

  def process_double(self, file_name):
    # Step 1: Load structures
    double_name = file_name.split('/')[-1].replace('.pdb', '') # _0001 should have been replaced befor thie step by some Shell script
    double = read_pdb(file_name)

    single_name_a, single_name_b = double_name.split('-')
    single_a = self.single_pdbs[single_name_a]
    single_b = self.single_pdbs[single_name_b]

    rcA = get_residue_count(single_a)
    rcB = get_residue_count(single_b)
    rc_double = get_residue_count(double)

    start_resi = int_floor(float(rcA)/2)
    rc_b_end = int_ceil(float(rcB)/2)
    end_resi = rc_double - rc_b_end

    #   The fusion_count is the number of residues we use to align double
    # to single_a. The high this number is, the more global our alignment
    # is, which causes bad disconnections in the chain. This is because in
    # Synth we're inevitably fusing atom positions from different doubles
    # into the same chain. Different doubles have their single components
    # stuck together using and interface, the participation of which
    # causes atom positions in the single components to differ from that
    # of the original single module.
    #
    #   When we fuse different doubles together, each double is cut at 25%
    # and 75% of their sequence in order to be as far way to interfaces
    # (0%, 50%, 100%) as possible.
    #
    #   The fusion alignment here is about aligning sub- sequent doubles
    # using a few residues before the 25% mark. The lower the fusion_count
    # is, the fewer resi- dues we use to align and the more local the
    # align- ment. However, if this number is too low the align- ment
    # could cause subsequent modules to overlap.
    #
    #   Through some experients I found that using 1/8 of the length of
    # single_b is a good balance between not causing discontinuities and
    # also not creating atom overlaps.
    fusion_count = int_ceil(float(rcB)/8)

    # Step 2: Move double to align with first single
    #
    #   This aligns double by superimposing double[0] with single_a
    #
    #   We only want to align to the second quardrant of single_a's atoms.
    # This is to be consistent with step 5
    self.align(
      double, 
      single_a, 
      match_count=start_resi
    )

    # Step 3: Get COM of the single_b as seen in the double
    #
    #    Only align the second quardrant of single_b in order to be
    # consistent with step 5
    comB = self.get_centre_of_mass(
      single_b, 
      mother=double, 
      child_resi_offset=rc_b_end - fusion_count,
      mother_resi_offset=rcA + rc_b_end - fusion_count,
      match_count=fusion_count
    )

    # Step 4: Get radius for collision checks later:
    #           1. Avg dist to com (gyradius aka RG)
    #           2. Max dist from CA to com
    #           3. Max dist from any heavy stom (not H) to COM
    radA = self.get_radii(single_a)

    # Step 5: Get transformation of double to the second single
    #
    #   Double is already aligned to first single so there is no need for
    # the first transformation
    #
    #   You can check this is true by varifying that
    # self.get_rot_trans(double, single_a) has identity rotation and zero
    # translation.
    #
    #   Only align the second quardrant of single_b in order to be
    # consistent with the Synth script, where doubles are fused together
    # by chopping the first and last quardrant of a double. This means the
    # second half of single_b is chopped off during fusion, while the first
    # quardrant of single_b participates in interfacing. Therefore we align
    # by uperimposing just the second quardrant
    rot, tran = self.get_rot_trans(
      double, 
      single_b, 
      moving_resi_offset=(rcA + rc_b_end - fusion_count),
      fixed_resi_offset=rc_b_end - fusion_count,
      match_count=fusion_count
    )

    # Step 6: Save the aligned molecules
    #
    #   Here the PDB format adds some slight floating point error. It is
    # really old and we should consider using mmCIF
    save_pdb(
      double, 
      self.aligned_pdb_dir + '/doubles/' + double_name + '.pdb'
    )

    data = OrderedDict([
      ('comB',  comB.tolist()),
      ('rot',   rot.tolist()),
      ('tran',  tran.tolist())
    ])

    entry = self.doubles_data.get(single_name_a, None)
    if entry == None:
      self.doubles_data[single_name_a] = {}
      entry = self.doubles_data.get(single_name_a, None)

    entry[single_name_b] = data

    single_data_a = self.singles_data.get(
      single_name_a,
      OrderedDict([
        ('link_count', 0),
        ('radii', radA)
      ])
    );
    single_data_a['link_count'] = single_data_a['link_count'] + 1;
    self.singles_data[single_name_a] = single_data_a;

  def dump_xdb(self):
    to_dump = OrderedDict([
      ('singles_data', self.singles_data),
      ('doubles_data', self.doubles_data)
      ])

    json.dump(to_dump,
      open(self.out_file, 'w'),
      separators=(',', ':'),
      ensure_ascii=False,
      indent=4)

  def run(self):
    # Center all single modules
    single_files = glob.glob(self.relaxed_pdbs_dir + '/singles/*.pdb')
    n_singles = len(single_files)
    for i in range(0, n_singles):
      print 'Centering single [{}/{}] {}' \
        .format(i+1, n_singles, single_files[i])
      self.process_single(single_files[i])

    # _0001 stands for relaxed PDBs
    doubleFiles = glob.glob(self.relaxed_pdbs_dir + '/doubles/*.pdb')
    nDoubles = len(doubleFiles)
    for i in range(0, nDoubles):
      print 'Aligning double [{}/{}] {}' \
        .format(i+1, nDoubles, doubleFiles[i])
      self.process_double(doubleFiles[i])

    print 'Total: {} singles, {} doubles'.format(n_singles, nDoubles)

    self.dump_xdb()

if __name__ =='__main__': 
  safe_exec(main)