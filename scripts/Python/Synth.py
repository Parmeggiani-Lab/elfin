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
import Bio.PDB
import Bio.SubsMat.MatrixInfo
import Bio.PDB.StructureBuilder
from ElfinUtils import *

class Synthesiser:
    def __init__(self, spec, pairsDir, singlesDir, crDir, showFusion=False):
        self.spec           = spec
        self.pairsDir       = pairsDir
        self.singlesDir     = singlesDir
        self.crDir          = crDir
        self.showFusion     = showFusion
        self.si             = Bio.PDB.Superimposer()

        # Parse and convert capping repeat indicies into a dictionary
        self.crIndexDict = {}
        for row in readCsv(self.crDir + '/repeat_indicies.csv', delim=' '):
            self.crIndexDict[row[0].split('.')[0].replace('DHR', 'D')] = [int(idx) for idx in row[1:]]

    def stripResidues(self, pdb):
        chain = getChain(pdb)
        residues = [r for r in chain.child_list]
        [chain.detach_child(r.id) for r in residues]
        return residues

    def resetResidueId(self):
        self.residueID = 1

    def getNextResidueId(self):
        rid = self.residueID
        self.residueID += 1
        return rid

    def capTerminus(self, primRes, capRes, crIndicies, nTerm=True):
        nCapRes = len(capRes)

        # See GenXDB for the division magic number
        # processPair() step 5
        if nTerm:
            replaceLen  = crIndicies[0] - 1
            alignResIdx = range(replaceLen + 1, (replaceLen+(nCapRes/8)) + 1)
        else:
            replaceLen  = len(capRes) - (crIndicies[3]-crIndicies[2]+1)
            alignResIdx = range(-replaceLen-(nCapRes/8) + 1, -replaceLen + 1)

        primAlignRes    = [primRes[i] for i in alignResIdx]
        capAlignRes     = [capRes[i] for i in alignResIdx]
        pauseCode()

        primAtoms       = [[a for a in r.child_list if a.name == 'CA'] for r in primAlignRes]
        capAtoms        = [[a for a in r.child_list if a.name == 'CA'] for r in capAlignRes]
        
        self.si.set_atoms(capAtoms, primAtoms)

        # Align and replace capRes onto primRes

    def projectModuleInstance(
        self, 
        graph,
        chains, 
        nodeA
    ):
        if not nodeA['trim'][1]:
            # This is an end-node and end-node atoms 
            # are already covered by their preceeding 
            # non-end node
            return chains

        print 'Processing node: id={}, name={}'.format(nodeA['id'], nodeA['name'])

        singleA = readPdb(self.singlesDir + '/' + nodeA['name'] + '.pdb')
        resiCountA = getResidueCount(singleA)

        ctermNodes = ctermNodes=[n for n in graph['nodes'] if n['id'] == nodeA['ctermNodeId']]
        nCtermNodes = len(ctermNodes)

        if nCtermNodes == 0:
            print ('Warning: \n'
                    '   Single-node chains still use transformation '
                    '   matricies calculated with pair fusion in mind. '
                    '   This might cause some inaccuracy when applied '
                    '   not on a trimmed pair but a single module.')
            chainId = graph['name'] + '_' + str(nodeA['id'])
            residues = self.stripResidues(singleA)

        elif nCtermNodes == 1:
            nodeB = ctermNodes[0]

            singleB = readPdb(self.singlesDir + '/' + nodeB['name'] + '.pdb')
            pair = readPdb(self.pairsDir + '/' + (nodeA['name'] + '-' + nodeB['name']) + '.pdb')

            resiCountB = getResidueCount(singleB)
            resiCountPair = getResidueCount(pair)

            chainId = graph['name'] + '_' + str(nodeA['id']) + '-' + str(nodeB['id'])
            residues = self.stripResidues(pair)

            # Compute trim residue indicies that minimise effect on 
            # atom position caused by pair interfaces
            if nodeA['trim'][0]:
                startResi = intFloor(float(resiCountA)/2)
            else:
                startResi = 0

                # We expect an untrimmed end to be capped
                if nodeA['cap'][0]:
                    crIndicies = self.crIndexDict[nodeA['name']]
                    cap = readPdb(self.crDir + '/' + nodeA['name'] + '_NI.pdb')
                    # self.capTerminus(
                    #     primRes=residues, 
                    #     capRes=self.stripResidues(cap), 
                    #     crIndicies=crIndicies, 
                    #     nTerm=True
                    # )
                    # pauseCode()
                else:
                    print 'Warning: untrimmed N-Terminus is not capped'
            
            if nodeB['trim'][1]:
                endResi = resiCountPair - intCeil(float(resiCountB)/2)
            else:
                endResi = resiCountPair

                # We expect an untrimmed end to be capped
                if nodeB['cap'][1]:
                    pauseCode()
                    crIndicies = self.crIndexDict[nodeB['name']]
                else:
                    print 'Warning: untrimmed C-Terminus is not capped'

            residues = residues[startResi:endResi]
        else:
            raise ValueError('nCtermNodes = ' + str(nCtermNodes))

        # Steal residues from the PDB and put in our chain
        if self.showFusion:
            currChain = Bio.PDB.Chain.Chain(chainId)
            chains.append(currChain)
        else:
            currChain = chains[0]

        for r in residues:
            r.id = (r.id[0], self.getNextResidueId(), r.id[2]) 
            r.transform(np.asarray(nodeA['rot']), nodeA['tran'])
            currChain.add(r)

    def makeChains(self, graph):
        chains = []

        if not self.showFusion:
            chains.append(Bio.PDB.Chain.Chain(graph['name']))

        for node in graph['nodes']:
            self.projectModuleInstance(
                graph=graph,
                chains=chains, 
                nodeA=node
            )

        return chains

    def run(self):
        model = Bio.PDB.Model.Model(0)
        si = Bio.PDB.Superimposer()

        if self.showFusion:
            print 'showFusion is on'

        self.resetResidueId()
        for graphIdx in xrange(0, len(self.spec)):
            print 'Processing graph: idx={}, name={}'.format(graphIdx, self.spec[graphIdx]['name'])
            map(lambda(c): model.add(c), self.makeChains(self.spec[graphIdx]))

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
    ap.add_argument('--cappingRepeatsDir', default='res/capping_repeats')
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
        print 'Warning: multi-chain feature is not well tested yet'

    model = Synthesiser(
        spec, 
        args.pairsDir,
        args.singlesDir,
        args.cappingRepeatsDir,
        args.showFusion
    ).run()

    if args.outFile == '':
        args.outFile = args.specFile
    args.outfile = '.'.join(args.outFile.split('.')[:-1])

    print 'Saving CIF...'
    saveCif(model, args.outFile + '.cif')

    # Todo: output coms for multiple chains
    #   Need to change csv format such that points
    # aren't connected between chains
    
    # coms = [map(lambda(el): el['tran'], spec[0]['nodes'])]
    # saveCsv(coms, args.outFile + '.csv')

if __name__ == '__main__':
    main()