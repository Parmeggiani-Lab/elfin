#!/usr/bin/env python
import glob, sys
from ElfinUtils import *
import Bio.PDB

### Sliding window based
if len(sys.argv) < 3:
    print './RMSDStat.py <solutionDir> <minimisedDir> <windowLen=300> <overlapRatio=0.50> <warnThreshold=5.0>'
    exit()

sPdbDir = sys.argv[1]
mPdbDir = sys.argv[2]
windowLen = 300
if len(sys.argv) >= 4:
	windowLen = int(sys.argv[3])
die(windowLen < 0, 'Invalid windowLen: must be greater than 0')

overlapRatio = 0.50
if len(sys.argv) >= 5:
	overlapRatio = float(sys.argv[4])
die(overlapRatio < 0 or overlapRatio > 1.0,
	'Invalid overlapRatio: must be between 0 and 1.0')

warnThreshold = 5.0
if len(sys.argv) >= 6:
	warnThreshold = float(sys.argv[5])
die(warnThreshold < 0,
	'Invalid warnThreshold: must be greater than 0')

overlap = int(overlapRatio * windowLen)
print 'Using window length {} and overlap ratio {} (={} CA atoms), warn at {}A'.format(
	windowLen, overlapRatio, overlap, warnThreshold)

mFiles = glob.glob(mPdbDir + '/*.pdb')

for i in range(0, len(mFiles)):
	mFile = mFiles[i]
	sFile = sPdbDir + mFile[mFile.rfind('/'):].replace('_0001.pdb', '.pdb')
		
	mPdb = readPdb(mFile)
	sPdb = readPdb(sFile)

	mAtoms = []
	sAtoms = []
	for a in mPdb.get_atoms():
		if a.id == 'CA':
			mAtoms.append(a)
	for a in sPdb.get_atoms():
		if a.id == 'CA':
			sAtoms.append(a)

	# Superimpose the two structures before comparing
	si = Bio.PDB.Superimposer()
	si.set_atoms(mAtoms, sAtoms)
	rot, tran = si.rotran

	sPdb.transform(rot, tran)

	mCAs = []
	sCAs = []
	for a in mPdb.get_atoms():
		if a.id == 'CA':
			mCAs.append(a.coord)
	for a in sPdb.get_atoms():
		if a.id == 'CA':
			sCAs.append(a.coord)

	mCAsLen =len(mCAs)

	if mCAsLen != len(sCAs):
		print('{}: Fatal! Number of CAs are different... Solution: {}, Minimised: {}'.format(
			sFile, len(sCAs), mCAsLen))
	else:
		stats = {
			'min': float('inf'),
			'avg': float('nan'),
			'max': -float('inf')
		}

		sumWinRmsd = 0.0
		startIndex = 0
		windowCount = 0
		while startIndex + windowLen - 1 < mCAsLen:
			sumD = 0.0
			for i in xrange(startIndex, startIndex+windowLen):
				sumD += np.linalg.norm(mCAs[i] - sCAs[i], 2)

			winRmsd = np.sqrt(sumD / windowLen)
			if winRmsd < stats['min']:
				stats['min'] = winRmsd
			if winRmsd > stats['max']:
				stats['max'] = winRmsd

			sumWinRmsd += winRmsd
			windowCount += 1

			startIndex += overlap

		stats['avg'] = sumWinRmsd / windowCount

		print '{:20} avg: {:10} min {:10} max {:10}'.format(
			sFile, stats['avg'], stats['min'], stats['max'])

		if stats['max'] > warnThreshold:
			print 'Warning: max window RMSD of {} exceeded!'.format(warnThreshold)
		# pauseCode()
