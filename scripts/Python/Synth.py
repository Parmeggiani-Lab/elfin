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

class Synthesiser:
    def __init__(self, spec, pairsDir, singlesDir, showFusion=False):
        self.spec = spec
        self.pairsDir = pairsDir
        self.singlesDir = singlesDir
        self.showFusion = showFusion

    def getChain(self, struct):
        return struct.child_list[0].child_dict['A']

    def stripResidues(self, chain):
        residues = [r for r in chain.child_list]
        [chain.detach_child(r.id) for r in residues]
        return residues

    def resetResidueId(self):
        self.residueID = 1

    def getNextResidueId(self):
        rid = self.residueID
        self.residueID += 1
        return rid

    def desendBranch(
        self, 
        graph, 
        chains, 
        prevNode, 
        currNodeId):
        prevNodeName = prevNode['name']
        singleA = readPdb(
            prevNodeName,
            self.singlesDir + '/' + prevNodeName + '.pdb'
        )

        currNodes = [n for n in graph['nodes'] if n['id'] == currNodeId]
        assert(len(currNodes) == 1)
        currNode = currNodes[0]
        currNodeName = currNode['name']
        singleB = readPdb(
            currNodeName,
            self.singlesDir + '/' + currNodeName + '.pdb'
        )

        pairName = prevNode['name'] + '-' + currNode['name']
        pair = readPdb(
            pairName,
            self.pairsDir + '/' + pairName + '.pdb'
        )

        resiCountA = getResidueCount(singleA)
        resiCountB = getResidueCount(singleB)
        resiCountPair = getResidueCount(pair)

        if prevNode['isStart']:
            # First pair: ignore leading trim
            startResi = 0
            endResi = resiCountPair - intCeil(float(resiCountB)/2)
            # Todo: if not connected to this chain's own end, then
            #       use capping
            # Todo: deal with atom overlaps if node is a hub
        elif len(currNode['ctermNodes']) == 0:
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

        pairResidues = self.stripResidues(self.getChain(pair))

        if self.showFusion:
            currChain = Bio.PDB.Chain.Chain(graph['name'] + str(currNodeId))
            chains.append(currChain)
        else:
            currChain = chains[0]

        for r in pairResidues[startResi:endResi]:
            r.id = (r.id[0], self.getNextResidueId(), r.id[2]) 
            r.transform(np.asarray(prevNode['rot']), prevNode['tran'])
            currChain.add(r)

        for cnId in currNode['ctermNodes']:
            self.desendBranch(graph, chains, currNode, cnId)

    def formChains(self, graph):
        chains = []

        if not self.showFusion:
            chains.append(Bio.PDB.Chain.Chain(graph['name'] + '0'))

        nodes = graph['nodes']

        startingNodes = [n for n in nodes if n['isStart']]
        assert(len(startingNodes) == 1)
        startingNode = startingNodes[0]

        for cnId in startingNode['ctermNodes']:
            self.desendBranch(graph, chains, startingNode, cnId)

        return chains

    def run(self):
        model = Bio.PDB.Model.Model(0)
        si = Bio.PDB.Superimposer()

        if self.showFusion:
            print 'showFusion is on'

        self.resetResidueId()
        for graphIdx in xrange(0, len(self.spec)):
            map(lambda(c): model.add(c), self.formChains(self.spec[graphIdx]))

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
    ap.add_argument('--showFusion', action='store_true')

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

    model = Synthesiser(
        spec, 
        args.pairsDir,
        args.singlesDir,
        args.showFusion
    ).run()

    if args.outFile == '':
        args.outFile = args.specFile
    args.outfile = '.'.join(args.outFile.split('.')[:-1])

    saveCif(model, args.outFile + '_synth.cif')

    coms = [map(lambda(el): el['tran'], spec[0]['nodes'])]
    saveCsv(coms, args.outFile + '_synth.csv')

if __name__ == '__main__':
    main()