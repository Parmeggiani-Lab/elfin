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
    def __init__(
        self, 
        spec, 
        pdbDir, 
        cappingsDir,
        metadataDir, 
        showFusion=False,
        disableCapping=False
    ):
        self.spec           = spec
        self.doublesDir     = pdbDir + "/doubles/"
        self.singlesDir     = pdbDir + "/singles/"
        self.crDir          = cappingsDir
        self.showFusion     = showFusion
        self.disableCapping = disableCapping

        self.si             = Bio.PDB.Superimposer()

        # Parse and convert capping repeat indicies into a dictionary
        self.cappingRepeatRIdDict = {}
        for row in readCsv(metadataDir + '/repeat_indicies.csv', delim=' '):
            self.cappingRepeatRIdDict[row[0].split('.')[0].replace('DHR', 'D')] = [int(idx) for idx in row[1:]]

    def resetResidueId(self):
        self.residueID = 1

    def getNextResidueId(self):
        rid = self.residueID
        self.residueID += 1
        return rid

    def getCapping(self, primRes, capRes, crRids, nTerm=True):
        if self.disableCapping:
            return []

        nCapRes = len(capRes)

        # We could use as many shared residues as possible for alignment but
        # that could create gaps (jumps) due to some of the residues in
        # primary being affected by interface. Therefore here we use one eigth
        # like GenXDB does (repeat index range / 4 is 1/8 of the module).
        if nTerm:
            matchStartIdx = [i for (i,el) in enumerate(capRes) if el.id[1] == crRids[0]][0]
            for cri in xrange(matchStartIdx, nCapRes):
                primIdx = cri - matchStartIdx
                if primRes[primIdx].resname != capRes[cri].resname: 
                        break
                maxMatchIdx = cri

            maxAlignLen = intFloor((maxMatchIdx - matchStartIdx) / 4)
            primAlignRes = primRes[:maxAlignLen]
            capAlignRes = capRes[matchStartIdx:matchStartIdx+maxAlignLen]
            realCapRes = capRes[:matchStartIdx]

            print 'N-Terminus capping align len: {}'.format(maxAlignLen)
        else:
            matchEndIdx = [i for (i,el) in enumerate(capRes) if el.id[1] == crRids[3]][0]
            for cri in reversed(xrange(0, matchEndIdx + 1)):
                primIdx = cri - matchEndIdx - 1
                if primRes[primIdx].resname != capRes[cri].resname: 
                        break
                minMatchIdx = cri

            maxAlignLen = intFloor((matchEndIdx - minMatchIdx) / 4)
            primAlignRes = primRes[-maxAlignLen:]
            capAlignRes = capRes[matchEndIdx-maxAlignLen+1:matchEndIdx+1]
            realCapRes = capRes[matchEndIdx+1:]

            print 'C-Terminus capping align len: {}'.format(maxAlignLen)

        for i in xrange(0, len(primAlignRes)):
            assert(primAlignRes[i].resname == capAlignRes[i].resname)

        primAtoms       = [al[0] for al in [[a for a in r.child_list if a.name == 'CA'] for r in primAlignRes]]
        capAtoms        = [al[0] for al in [[a for a in r.child_list if a.name == 'CA'] for r in capAlignRes]]
        
        self.si.set_atoms(primAtoms, capAtoms)
        rot, tran = self.si.rotran

        [r.transform(rot, tran) for r in realCapRes]
        
        return realCapRes

    def projectModuleInstance(
        self, 
        graph,
        chains, 
        nodeA
    ):
        print 'Processing node: id={}, name={}'.format(nodeA['id'], nodeA['name'])

        if not nodeA['trim'][1]:
            # This is an end-node and end-node atoms are already covered by
            # their preceeding non-end node
            return chains

        singleA = readPdb(self.singlesDir + '/' + nodeA['name'] + '.pdb')
        resiCountA = getResidueCount(singleA)

        ctermNodes = ctermNodes=[n for n in graph['nodes'] if n['id'] == nodeA['ctermNodeId']]
        nCtermNodes = len(ctermNodes)

        if nCtermNodes == 0:
            # end-nodes don't reach here - this only fires in cases where
            # there's only one single node in the entire chain
            print ('Warning: \n'
                    '   Single-node chains still use transformation '
                    '   matricies calculated with double fusion in mind. '
                    '   This might cause some inaccuracy when applied '
                    '   not on a trimmed double but a single module.')
            chainId = graph['name'] + '_' + str(nodeA['id'])
            residues = stripResidues(singleA)

        elif nCtermNodes == 1:
            nodeB = ctermNodes[0]

            singleB = readPdb(self.singlesDir + '/' + nodeB['name'] + '.pdb')
            double = readPdb(self.doublesDir + '/' + (nodeA['name'] + '-' + nodeB['name']) + '.pdb')

            resiCountB = getResidueCount(singleB)
            resiCountDouble = getResidueCount(double)

            chainId = graph['name'] + '_' + str(nodeA['id']) + '-' + str(nodeB['id'])
            residues = stripResidues(double)

            # Compute trim residue indicies that minimise effect on 
            # atom position caused by double interfaces
            prefixResidues = []
            suffixResidues = []
            if nodeA['trim'][0]:
                startResi = intFloor(float(resiCountA)/2)
            else:
                startResi = 0

                # We expect an untrimmed end to be capped
                if nodeA['cap'][0]:
                    capName = nodeA['name'].split('_')[0]
                    capAndRepeat = readPdb(self.crDir + '/' + capName + '_NI.pdb')
                    prefixResidues = self.getCapping(
                        primRes=residues, 
                        capRes=stripResidues(capAndRepeat), 
                        crRids=self.cappingRepeatRIdDict[capName], 
                        nTerm=True
                    )
                else:
                    print 'Warning: untrimmed N-Terminus is not capped'
            
            if nodeB['trim'][1]:
                endResi = resiCountDouble - intCeil(float(resiCountB)/2)
            else:
                endResi = resiCountDouble

                # We expect an untrimmed end to be capped
                if nodeB['cap'][1]:
                    capName = nodeB['name'].split('_')[-1]
                    capAndRepeat = readPdb(self.crDir + '/' + capName + '_IC.pdb')
                    suffixResidues = self.getCapping(
                        primRes=residues, 
                        capRes=stripResidues(capAndRepeat), 
                        crRids=self.cappingRepeatRIdDict[capName], 
                        nTerm=False
                    )
                else:
                    print 'Warning: untrimmed C-Terminus is not capped'

            residues = prefixResidues + residues[startResi:endResi] + suffixResidues
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
            map(lambda c: model.add(c), self.makeChains(self.spec[graphIdx]))

        sb = Bio.PDB.StructureBuilder.StructureBuilder()
        sb.init_structure('0')
        structure = sb.get_structure()
        structure.add(model)
        return structure


def main():
    ap = argparse.ArgumentParser(description='Generate atom model in CIF format using output from Elfin core');
    ap.add_argument('specFile')
    ap.add_argument('--outFile', default='')
    ap.add_argument('--pdbDir', default='./resources/pdb_aligned/')
    ap.add_argument('--cappingsDir', default='./resources/pdb_relaxed/cappings')
    ap.add_argument('--metadataDir', default='./resources/metadata/')
    ap.add_argument('--showFusion', action='store_true')
    ap.add_argument('--disableCapping', action='store_true')
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
        args.pdbDir,
        args.cappingsDir,
        args.metadataDir,
        args.showFusion,
        args.disableCapping
    ).run()

    if args.outFile == '':
        args.outFile = args.specFile
    args.outfile = '.'.join(args.outFile.split('.')[:-1])

    print 'Saving...'
    saveCif(model, file)

    # Todo: output coms for multiple chains
    #   Need to change csv format such that points
    # aren't connected between chains
    
    # coms = [map(lambda el: el['tran'], spec[0]['nodes'])]
    # saveCsv(coms, args.outFile + '.csv')

if __name__ == '__main__':
    safeExec(main)