#!/usr/bin/env python
from ElfinUtils import *

def main():
	(avgD, minD, maxD) = getXDBStat(readJSON('resources/xDB.json'))
	print 'Distances avg: {}, min: {}, max: {}'.format(avgD, minD, maxD)

if __name__ =='__main__': safeExec(main)