#!/usr/bin/env python3

from utilities import *

def main():
	(avg_d, min_d, max_d) = com_dist_info(read_json('resources/xDB.json'))
	print('Distances avg: {}, min: {}, max: {}'.format(avgD, minD, maxD))

if __name__ =='__main__': 
  safe_exec(main)