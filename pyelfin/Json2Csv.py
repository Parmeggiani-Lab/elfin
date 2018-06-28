#!/usr/bin/env python
import glob, json, sys
from ElfinUtils import *

### This is a tool for converting spec JSON files into pure CoMs in CSV form ###

if len(sys.argv) < 2:
    print './Json2Csv.py <sourceDir> <outDir=sourceDir>'
    exit()

sourceDir = sys.argv[1]
outDir = sourceDir
if len(sys.argv) >= 3:
	outDir = sys.argv[2]

jsonFiles = glob.glob(sourceDir + '/*.json')

for i in range(0, len(jsonFiles)):
	inFile = jsonFiles[i]
	outFile = outDir + inFile[inFile.rfind('/'):].replace('.json', '.csv')

	inSpec = readJSON(inFile)

	outStr = ''
	for com in inSpec['coms']:
		outStr += '{} {} {}\n'.format(com[0], com[1], com[2])
	
	with open(outFile, 'w') as of:
		of.write(outStr)

print 'Done'