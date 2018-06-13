#!/usr/bin/env python

# 
# In this preprocessing stage, all doubles are chain-merged as they come in
# separate chains, which denote physically separate entities.
# 
# Then, unwanted atoms are deleted from each single and double.
#
# For "compound" doubles (those with junctions), loop interfaces must be
# replaced with those found in their "simple" counterparts (those without
# junctions).
#
# This needs to be done because we have experimental data that confirm
# interface atom positions in simple doubles but not in compound doubles. If
# this is not done, the relaxation of database PDBs will result in sever
# bending and dislocation in compound doubles.
#

import argparse, sys
from ElfinUtils import *

def mergeChainsAndCleanse(pdb):
	newChain = Bio.PDB.Chain.Chain('A')
	rid = 1
	for r in stripResidues(pdb):
		r.id = (r.id[0], rid, r.id[2])
		badAtoms = []
		for a in r:
			if a.name == '1H' or a.name == '2H' or a.name == '3H' or a.name == 'OXT':
				badAtoms.append(a.name)
		for ba in badAtoms:
			r.detach_child(ba)
		newChain.add(r)
		rid += 1

	# Remove old chains from model
	for model in pdb:
		ocIds = []
		for chain in model:
			ocIds.append(chain.id)
		for ocId in ocIds:
			model.detach_child(ocId)

	model.add(newChain)

def main():
	ap = argparse.ArgumentParser(description='Template Python script');
	ap.add_argument('input')
	ap.add_argument('--outputDir', default='./resources/pdb_prepped/')
	args = ap.parse_args()

	if args.input is None or args.outputDir is None:
		ap.print_help()
		sys.exit(1)

	# Extract name
	doubleFile = args.input
	doubleName = doubleFile[doubleFile.rfind('/')+1:].replace('.pdb', '')
	underscores = [doubleName.find('_'), doubleName.rfind('_')]

	if underscores[0] == -1:
		print 'Input {} is a simple double and does not need loop replacement'.format(args.input)
		double = readPdb(doubleFile)
		mergeChainsAndCleanse(double)
		savePdb(double, args.outputDir + '/' + doubleName + '.pdb')
		exit(0)

	dashIdx = doubleName.rfind('-')
	sdoubleFirst = dashIdx < underscores[0] and dashIdx < underscores[1]

	doubleNameHalves = [doubleName[:dashIdx], doubleName[dashIdx+1:]]
	sdoubleName = ''
	sdoubleName += doubleNameHalves[0][doubleNameHalves[0].rfind('_')+1:] \
				if doubleNameHalves[0].rfind('_') != -1 else doubleNameHalves[0]
	sdoubleName += '-'
	sdoubleName += doubleNameHalves[1][:doubleNameHalves[1].find('_')] \
				if doubleNameHalves[1].find('_') != -1 else doubleNameHalves[1]
	sdoubleFile = doubleFile[:doubleFile.rfind('/')+1] + sdoubleName + '.pdb'

	# Load PDBs
	double = readPdb(doubleFile)
	doubleChains = double.child_list[0].child_list
	assert(len(doubleChains) == 2)

	sdouble = readPdb(sdoubleFile) #sdouble is the simple double
	sdoubleChains = sdouble.child_list[0].child_list
	assert(len(sdoubleChains) == 2)

	# Get residue counts
	sdoubleRCount = getResidueCount(sdouble)

	sdoubleStartIdx = int(np.ceil(sdoubleRCount*0.375))+1
	sdoubleEndIdx = int(np.floor(sdoubleRCount*0.625))-1
	sdoubleEndOffset = sdoubleEndIdx - len(sdoubleChains[0].child_list)

	# Find simple double middle residues
	sdoubleChainLens = [len(c.child_list) for c in sdoubleChains]
	sdoubleMidRes = sdoubleChains[0].child_list[sdoubleStartIdx:] + \
		sdoubleChains[1].child_list[:sdoubleEndOffset]
	sdoubleAtoms = [a for r in sdoubleMidRes for a in r.child_list if a.name == 'CA']

	# Find first half of the residues
	doubleChainLens = [len(c.child_list) for c in doubleChains]
	doubleStartOffset = doubleChainLens[0]-(sdoubleChainLens[0]-sdoubleStartIdx)
	doubleMidRes = \
		(doubleChains[0].child_list[sdoubleStartIdx:] if sdoubleFirst else \
		doubleChains[0].child_list[doubleStartOffset:]) + \
		doubleChains[1].child_list[:sdoubleEndOffset]
	doubleAtoms = [a for r in doubleMidRes for a in r.child_list if a.name == 'CA']

	# Superimpose double onto sdouble
	si = Bio.PDB.Superimposer()
	si.set_atoms(sdoubleAtoms, doubleAtoms)
	rot, tran = si.rotran	
	double.transform(rot, tran)

	# Merge chains and remove bad atoms
	mergeChainsAndCleanse(double)
	mergeChainsAndCleanse(sdouble)

	# Replace double residues where it should be sdouble residues
	doubleResidues = getChain(double).child_list
	sdoubleResidues = stripResidues(sdouble)
	for rIdx in xrange(sdoubleStartIdx, sdoubleEndIdx):
		offsetRIdx = rIdx + (0 if sdoubleFirst else doubleChainLens[0]-sdoubleChainLens[0])
		oldRId = doubleResidues[offsetRIdx].id
		doubleResidues[offsetRIdx] = sdoubleResidues[rIdx]
		doubleResidues[offsetRIdx].id = oldRId

	savePdb(double, args.outputDir + '/' + doubleName + '.pdb')

if __name__ == '__main__':
	main()