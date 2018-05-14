#!/usr/bin/env python
import glob, sys
from ElfinUtils import *

### Rosetta overall score based... deprecated in favour of windowed RMSD (RMSDStat.py)
if(len(sys.argv) < 2):
    print './RMSDStat.py <scoreDir>'
    exit()

scoreDir = sys.argv[1]

files = glob.glob(scoreDir + '/*_comp.sc')
nFiles = len(files)
rmsds = []

for i in range(0, nFiles):
    with open(files[i], 'r') as file:
		line = file.read().split('\n')[-2]
		rmsdStr = line.split(' ')[-2]
		print '{} RMSD: {}'.format(files[i], rmsdStr)
		rmsds.append(float(rmsdStr))

maxRmsd = max(rmsds)
print 'Average: {}, Min: {}, Max: {}'.format(sum(rmsds)/nFiles, min(rmsds), maxRmsd)

if(maxRmsd > 5.0):
	print 'WARNING: One or more molecules exceed 5A RMSD!'