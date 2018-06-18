#!/usr/bin/env python

import argparse, glob, itertools
from ElfinUtils import *

def main():
  ap = argparse.ArgumentParser(description='Generate all combinations of hubs with other modules')
  ap.add_argument('--hubDir', default='./resources/pdb_raw/hubs/')
  ap.add_argument('--doublesDir', default='./resources/pdb_raw/doubles/')
  ap.add_argument('--outputDir', default='./resources/pdb_preppd/hubs/')
  ap.add_argument('--xdbPath', default='resources/xDB.json')
  args = ap.parse_args()

  doubleFiles = glob.glob(args.doublesDir + '/*.pdb')
  doubleNames = [os.path.basename(df).replace('.pdb', '') for df in doubleFiles]
  
  hub_info=readCsv('./resources/pdb_raw/hubs/hub_info.csv')
  for hb in glob.glob(args.hubDir + '/*.pdb'):
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

    print('Generating for hub {}'.format(hbName))

    # Compute all possible combinations

    # When D-type hubs come into the DB, the following might need to be
    # implemented:
    #   1. Indication of n_free and c_free PER chain (instead of universal)
    #   2. ?
    ntermDoubles = [None]
    ctermDoubles = [None]
    hbSingleName = info[5]
    if info[3].lower() == 'n_free':
      ntermDoubles += [dn for dn in doubleNames if dn.endswith('-' + hbSingleName)]
      print('Compatible N-term doubles: {}'.format(ntermDoubles))
    else:
      print('N-term not free')

    if info[4].lower() == 'c_free':
      ctermDoubles += [dn for dn in doubleNames if dn.startswith(hbSingleName + '-')]
      print('Compatible C-term doubles: {}'.format(ctermDoubles))
    else:
      print('C-term not free')

    nChains = int(info[6])

    symmetricity = info[2]
    if symmetricity == 'symmetric':
      setConfigs = list(itertools.product(ntermDoubles, ctermDoubles))
    elif symmetricity == 'asymmetric':
      ntermProductArgs = [ntermDoubles] * nChains
      ntermConfigs = list(itertools.product(*ntermProductArgs))
      ctermProductArgs = [ctermDoubles] * nChains
      ctermConfigs = list(itertools.product(*ctermProductArgs))
      setConfigs = list(itertools.product(ntermConfigs, ctermConfigs))
    else:
      raise ValueError('Invalid symmetricity field: \'{}\''.format(symmetricity))

    for sc in setConfigs:
      print(sc)
    print('Configs: {}'.format(len(setConfigs)))
    pauseCode()


if __name__ == '__main__':
  safeExec(main)