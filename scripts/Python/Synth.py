#!/usr/bin/env python

#
# This script creates the atom model of an Elfin design solution
#   Input: a JSON file that describes the connectivity of a sol-
#       lution. 
#   Output: a CIF file define the positions of each atom of the 
#       input solution
#
# Version: v2
#
#

import glob
import argparse
import Bio.SubsMat.MatrixInfo
import Bio.PDB.StructureBuilder
from ElfinUtils import *

def getChain(struct):
    return struct.child_list[0].child_dict['A']

def stripResidues(chain):
    residues = [r for r in chain.child_list]
    [chain.detach_child(r.id) for r in residues]
    return residues

def synthesise(xdb, spec, pairsDir, singlesDir, forceSplitChains=False):
    model = Bio.PDB.Model.Model(0)
    si = Bio.PDB.Superimposer()

    if forceSplitChains:
        print 'forceSplitChains is on'

    for graphIdx in xrange(0, len(spec)):
        graph = spec[graphIdx]
        currChain = Bio.PDB.Chain.Chain(str(graphIdx))

        nodes = graph['nodes']
        chainLenDigits = len(str(len(nodes)))

        residueUid = 1
        # firstNodeStruct = readPdb(
        #     nodes[0]['name'],
        #     singlesDir + '/' + nodes[0]['name'] + '.pdb'
        # )
        # for r in stripResidues(getChain(firstNodeStruct)):
        #     r.id = (r.id[0], residueUid, r.id[2]) 
        #     currChain.add(r)
        #     residueUid += 1

        for nodeIdx in xrange(1, len(nodes)):
            # Todo: use edges instead of looping thorugh nodes

            prevNode = nodes[nodeIdx-1]
            currNode = nodes[nodeIdx]
            rel = xdb['pairsData'][prevNode['name']][currNode['name']]

            # currChain.transform(np.asarray(rel['rot']), rel['tran'])
            # currNodeStruct = readPdb(
            #     currNode['name'],
            #     singlesDir + '/' + currNode['name'] + '.pdb'
            # )
            # for r in stripResidues(getChain(currNodeStruct)):
            #     r.id = (r.id[0], residueUid, r.id[2]) 
            #     currChain.add(r)
            #     residueUid += 1

            # if forceSplitChains:
            #     model.transform(np.asarray(rel['rot']), rel['tran'])
            #     model.add(currChain.copy())
            #     currChain._id = str(nodeIdx)
            #     stripResidues(currChain)
            
            prevNodeName = prevNode['name']
            singleA = readPdb(
                prevNodeName,
                singlesDir + '/' + prevNodeName + '.pdb'
            )
            currNodeName = currNode['name']
            singleB = readPdb(
                currNodeName,
                singlesDir + '/' + currNodeName + '.pdb'
            )
            pairName = prevNode['name'] + '-' + currNode['name']
            pair = readPdb(
                pairName,
                pairsDir + '/' + pairName + '.pdb'
            )

            resiCountA = getResidueCount(singleA)
            resiCountB = getResidueCount(singleB)
            resiCountPair = getResidueCount(pair)

            if nodeIdx == 1:
                # First pair: ignore leading trim
                startResi = 0
                endResi = resiCountPair - intCeil(float(resiCountB)/2)
                # Todo: if not connected to this chain's own end, then
                #       use capping
            elif nodeIdx == len(nodes) - 1:
                # Last pair: ignore trailing trim
                startResi = intFloor(float(resiCountA)/2)
                endResi = resiCountPair
                # Todo: if not connected to this chain's own begining,
                #       then use capping
            else:
                # Trim half of singleA's residues and extend
                # to half of singleB's residues
                startResi = intFloor(float(resiCountA)/2)
                endResi = resiCountPair - intCeil(float(resiCountB)/2)

            pairResidues = stripResidues(getChain(pair))

            for r in pairResidues[startResi:endResi]:
                r.id = (r.id[0], residueUid, r.id[2]) 
                currChain.add(r)
                residueUid += 1
            
            currChain.transform(np.asarray(rel['rot']), rel['tran'])

        if not forceSplitChains:
            model.add(currChain)

    sb = Bio.PDB.StructureBuilder.StructureBuilder()
    sb.init_structure('0')
    structure = sb.get_structure()
    structure.add(model)
    return structure

def main():
    ap = argparse.ArgumentParser(description='Generate atom model in CIF format using output from Elfin core');
    ap.add_argument('specFile')
    ap.add_argument('--outFile', default='')
    ap.add_argument('--singlesDir', default='res/aligned/single/')
    ap.add_argument('--pairsDir', default='res/aligned/pair/')
    ap.add_argument('--xdbPath', default='res/xDB.json')
    ap.add_argument('--forceSplitChains', action='store_true')

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)
        
    args = ap.parse_args()

    specExt = args.specFile[args.specFile.rfind('.'):]

    if specExt == '.json':
        spec = readJSON(args.specFile)
    else:
        print 'Unknown spec file type: {}'.format(specExt)
        exit()

    if len(spec) > 1:
        print 'Not yet implemented: fill atoms for multiple chains'
        exit()

    xDB = readJSON(args.xdbPath)
    model = synthesise(
        xDB, 
        spec, 
        args.pairsDir,
        args.singlesDir,
        args.forceSplitChains
    )

    if args.outFile == '':
        args.outFile = args.specFile
    args.outfile = '.'.join(args.outFile.split('.')[:-1])

    saveCif(model, args.outFile + '_synth.cif')

    coms = [map(lambda(el): el['tran'], spec[0]['nodes'])]
    saveCsv(coms, args.outFile + '_synth.csv')

if __name__ == '__main__':
    main()