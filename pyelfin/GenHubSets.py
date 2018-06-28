#!/usr/bin/env python

import argparse, glob, itertools
from ElfinUtils import *

def main():
  ap = argparse.ArgumentParser(description='Generate all combinations of hubs with other modules')
  ap.add_argument('--pdbDir', default='./resources/pdb_relaxed/')
  ap.add_argument('--outputDir', default='./resources/pdb_preppd/hubs/')
  ap.add_argument('--xdbPath', default='resources/xDB.json')
  args = ap.parse_args()

  doubleFiles = glob.glob(args.pdbDir + '/doubles/*.pdb')
  doubleNames = [os.path.basename(f).replace('.pdb', '') for f in doubleFiles]
  
  totalSetConfigs = 0
  hub_info=readCsv('./resources/metadata/hub_info.csv')
  for hb in glob.glob(args.pdbDir + '/hubs/*.pdb'):
    hbBasename = os.path.basename(hb)
    hbName = hbBasename.replace('.pdb', '')

    # Check that the current hub has corresponding info row
    try:
      info = [line for line in hub_info if line[0] == hbBasename][0]
    except Exception as e:
      print("Error: cannot find hub inside hub info data")
      print('Basename: {}'.format(hbBasename))
      for line in hub_info:
        print('Line=\'{}\', found={}'.format(line, line[0] == hbBasename))
      raise e

    print('Generating hub set for {}'.format(hbName))

    # Compute all possible combinations

    # When D-type hubs come into the DB, the following might need to be
    # implemented:
    #   1. Indication of n_free and c_free PER chain (instead of universal)
    #   2. ?
    ntermModules = [None]
    ctermModules = [None]
    hbSingleName = info[5]
    if info[3].lower() == 'n_free':
      ntermModules = list(set([d.split('-')[0] for d in doubleNames if ('-' + hbSingleName) in d]))
      print('Compatible N-term modules: {}'.format(ntermModules))
    else:
      print('N-term not free')
      ntermModules = [None]

    if info[4].lower() == 'c_free':
      ctermModules += list(set([d.split('-')[1] for d in doubleNames if (hbSingleName + '-') in d]))
      print('Compatible C-term modules: {}'.format(ctermModules))
    else:
      print('C-term not free')
      ctermModules = [None]

    nChains = int(info[6])

    symmetricity = info[2]
    if symmetricity == 'symmetric':
      print('Hub is symmetric')
      setConfigs = list(itertools.product([(m,)*nChains for m in ntermModules], [(m,)*nChains for m in ctermModules]))
    elif symmetricity == 'asymmetric':
      print('Hub is asymmetric')
      ntermProductArgs = [ntermModules] * nChains
      ntermConfigs = list(itertools.product(*ntermProductArgs))
      ctermProductArgs = [ctermModules] * nChains
      ctermConfigs = list(itertools.product(*ctermProductArgs))
      setConfigs = list(itertools.product(ntermConfigs, ctermConfigs))
    else:
      raise ValueError('Invalid symmetricity field: \'{}\''.format(symmetricity))

    # Do the maths
      
    # from ElfinUtils import *
    # hub = readPdb('/mnt/c/Users/Akaoni/Desktop/ElfinWork/elfin/resources/pdb_raw/hubs/D4_C3_02_0001.pdb')
    # double = readPdb('/mnt/c/Users/Akaoni/Desktop/ElfinWork/elfin/resources/pdb_relaxed/doubles/D4-D4.pdb')

    # hubA = getChains(hub)[0]
    # dbChain = getChains(double)[0]
    # dbOffset = int(len(dbChain.child_list)/2)

    # hubACAs = [a for r in hubA.child_list[:3] for a in r.child_list if a.name == 'CA']
    # doubleBCAs = [a for r in dbChain.child_list[dbOffset:dbOffset+3] for a in r.child_list if a.name == 'CA']

    # si = Bio.PDB.Superimposer()
    # si.set_atoms(hubACAs, doubleBCAs)
    # double.transform(*si.rotran)

    # hubARs = stripResidues(hub, chainIds=['A'])
    # doubleRs = stripResidues(double, chainIds=['A'])[:dbOffset]
    # ARs = doubleRs + hubARs

    # rid = 1
    # for r in ARs:
    #   r.id = (r.id[0], rid, r.id[2])
    #   rid += 1
    #   hubA.add(r)

    # saveCif(hub, '/mnt/c/Users/Akaoni/Desktop/ElfinWork/elfin/resources/pdb_raw/hubs/test.cif')

    totalSetConfigs += len(setConfigs);
    for sc in setConfigs:
      print(sc)
    print('Configs: {}'.format(len(setConfigs)))
  print('Total configs: {}'.format(totalSetConfigs))

if __name__ == '__main__':
  safeExec(main)