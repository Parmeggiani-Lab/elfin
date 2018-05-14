#!/usr/bin/env python

# In this preprocessing stage, all pairs are chain-merged 
# as they come in separate chains, which denote physically
# separate entities.
# 
# Then, unwanted atoms are deleted from each single and pair.
#
# For "compound" pairs (those with junctions), loop interfaces 
# are replaced with those found in their "simple" pair counter-
# parts (those withou junctions).
#
# This needs to be done because we are only sure about inter-
# face atom positions in simple pairs but not in compound 
# pairs. If this is not done, the relaxation of database PDBs 
# will result in sever bending and dislocation in compound pairs.
#

import argparse, sys
from ElfinUtils import *

def mergeChainsAndCleanse(pdb):
	for model in pdb:
		newChain = Bio.PDB.Chain.Chain('A')
		rid = 1
		for chain in model:
			for r in chain:
				badAtoms = []
				for a in r:
					if a.name == '1H' or a.name == '2H' or a.name == '3H' or a.name == 'OXT':
						badAtoms.append(a.name)
				for ba in badAtoms:
					r.detach_child(ba)


				r.id = (r.id[0], rid, r.id[2])
				newChain.add(r)
				rid += 1


		# Remove old chains from model
		ocIds = []
		for chain in model:
			ocIds.append(chain.id)

		for ocId in ocIds:
			model.detach_child(ocId)

		model.add(newChain)

def main():
	ap = argparse.ArgumentParser(description='Template Python script');
	ap.add_argument('input')
	ap.add_argument('--outputDir', default='./res/preprocessed/')
	args = ap.parse_args()

	if args.input is None or args.outputDir is None:
		ap.print_help()
		sys.exit(1)

	# Extract name
	pairFile = args.input
	pairName = pairFile[pairFile.rfind('/')+1:].replace('.pdb', '')
	underscores = [pairName.find('_'), pairName.rfind('_')]

	print 'Working on {}'.format(args.input)
	if underscores[0] == -1:
		print 'Input {} is a simple pair and does not need loop replacement'.format(args.input)
		pair = readPdb('pair', pairFile)
		mergeChainsAndCleanse(pair)
		savePdb(pair, args.outputDir + '/' + pairName + '.pdb')
		exit(0)

	dashIdx = pairName.rfind('-')
	spairFirst = dashIdx < underscores[0] and dashIdx < underscores[1]

	pairNameHalves = [pairName[:dashIdx], pairName[dashIdx+1:]]
	spairName = ''
	spairName += pairNameHalves[0][pairNameHalves[0].rfind('_')+1:] \
				if pairNameHalves[0].rfind('_') != -1 else pairNameHalves[0]
	spairName += '-'
	spairName += pairNameHalves[1][:pairNameHalves[1].find('_')] \
				if pairNameHalves[1].find('_') != -1 else pairNameHalves[1]
	spairFile = pairFile[:pairFile.rfind('/')+1] + spairName + '.pdb'

	# Load PDBs
	pair = readPdb('pair', pairFile)
	pairChains = pair.child_list[0].child_list
	assert(len(pairChains) == 2)

	spair = readPdb('spair', spairFile) #spair is the simple pair
	spairChains = spair.child_list[0].child_list
	assert(len(spairChains) == 2)

	# Get residue counts
	pairRCount = getResidueCount(pair)
	spairRCount = getResidueCount(spair)

	spairStartIdx = int(np.ceil(spairRCount*0.375))+1
	spairEndIdx = int(np.floor(spairRCount*0.625))-1
	spairEndOffset = spairEndIdx - len(spairChains[0].child_list)

	# Find simple pair middle residues
	spairChainLens = [len(c.child_list) for c in spairChains]
	spairMidRes = spairChains[0].child_list[spairStartIdx:] + \
		spairChains[1].child_list[:spairEndOffset]
	spairAtoms = [a for r in spairMidRes for a in r.child_list if a.name == 'CA']

	# Find first half of the residues
	pairChainLens = [len(c.child_list) for c in pairChains]
	pairStartOffset = pairChainLens[0]-(spairChainLens[0]-spairStartIdx)
	pairMidRes = \
		(pairChains[0].child_list[spairStartIdx:] if spairFirst else \
		pairChains[0].child_list[pairStartOffset:]) + \
		pairChains[1].child_list[:spairEndOffset]
	pairAtoms = [a for r in pairMidRes for a in r.child_list if a.name == 'CA']

	# Superimpose pair onto spair
	si = Bio.PDB.Superimposer()
	si.set_atoms(spairAtoms, pairAtoms)
	rot, tran = si.rotran	
	pair.transform(rot, tran)

	# Merge chains and remove bad atoms
	mergeChainsAndCleanse(pair)
	mergeChainsAndCleanse(spair)

	# Replace pair residues where it should be spair residues
	pairNewChain = pair.child_list[0].child_list[0]
	spairNewChain = spair.child_list[0].child_list[0]
	for rIdx in xrange(spairStartIdx, spairEndIdx):
		offsetRIdx = rIdx + (0 if spairFirst else pairChainLens[0]-spairChainLens[0])
		oldRId = pairNewChain.child_list[offsetRIdx].id
		pairNewChain.child_list[offsetRIdx] = spairNewChain.child_list[rIdx]
		pairNewChain.child_list[offsetRIdx].id = oldRId

	savePdb(pair, args.outputDir + '/' + pairName + '.pdb')

if __name__ == '__main__':
	main()