#!/usr/bin/env python3

import argparse, sys
from collections import OrderedDict
from utilities import *

def main():
  ap = argparse.ArgumentParser(description='Converts hub info metadata from csv to json.')
  ap.add_argument('input') # Absence of dash denotes mandatory argument
  args = ap.parse_args()

  csv_data = read_csv(args.input)
  new_data = OrderedDict({})
  for i in range(int(len(csv_data)/2)):
    row = csv_data[int(2*i)]
    chains = csv_data[int(2*i+1)]
    component_info = OrderedDict({})
    for i in range(int(len(chains)/2)):
      component_info[chains[int(2*i)]] = { 'single_name': chains[int(2*i+1)], 'c_free': row[4]=='C_free', 'n_free': row[3]=='N_free' }
    new_data[row[0].replace('.pdb', '')] = OrderedDict({ 'oligomer_type': row[1], 'symmetric': row[2] == 'symmetric', 'component_info': component_info})
  json.dump(new_data, open(args.input.replace('.csv', '') + '.json', 'w'), separators=(',',':'), ensure_ascii=False, indent=4)

if __name__ == '__main__':
  safe_exec(main)