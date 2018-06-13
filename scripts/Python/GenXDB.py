#!/usr/bin/env python

import Bio.PDB
import ElfinUtils
import glob
import numpy as np
import codecs, json
import argparse
from collections import OrderedDict

def main():
    ap = argparse.ArgumentParser(description='Generates the xDB database from preprocessed single and pair modules.');
    ap.add_argument('--outputAlignedDir', default='./res/aligned_modules/')
    ap.add_argument('--pairDir', default='./res/relaxed_modules/pair/')
    ap.add_argument('--singleDir', default='./res/relaxed_modules/single/')
    ap.add_argument('--output', default='./res/xDB.json')
    ap.add_argument('--genFakeHubs', action='store_true')
    args = ap.parse_args()

    xdbg = XDBGenrator(
        args.pairDir, 
        args.singleDir, 
        args.outputAlignedDir, 
        args.output, 
        args.genFakeHubs)
    ElfinUtils.safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                outputAlignedDir,
                outFile,
                genFakeHubs):
        self.pairDir        = pairDir
        self.singleDir      = singleDir
        ElfinUtils.mkdir(outputAlignedDir)
        ElfinUtils.mkdir(outputAlignedDir     + '/pair/')
        ElfinUtils.mkdir(outputAlignedDir     + '/single/')
        self.outputAlignedDir  = outputAlignedDir
        self.outFile        = outFile
        self.genFakeHubs    = genFakeHubs
        self.si             = Bio.PDB.Superimposer()
        self.pairsData      = {}
        self.singlesData    = {}
        self.singlePDBs     = {}

    def getCOM(
            self, 
            child, 
            mother=None, 
            childResiOffset=0, 
            motherResiOffset=0,
            matchCount=-1
        ):
        CAs = []
        for a in child.get_atoms():
            if(a.name == 'CA'):
                CAs.append(a.get_coord().astype('float64'))
        com = np.mean(CAs, axis=0)

        if mother is not None:
            # This is for finding COM of a single inside a pair
            _, tran = self.getRotTrans(
                child, 
                mother, 
                movingResiOffset=childResiOffset, 
                fixedResiOffset=motherResiOffset,
                matchCount=matchCount
            )

            com += tran

        return com

    def moveToOrigin(self, pdb):
        com = self.getCOM(pdb)

        # No rotation - just move to centre
        pdb.transform([[1,0,0],[0,1,0],[0,0,1]], -com)

    def align(
            self, 
            moving, 
            fixed, 
            movingResiOffset=0, 
            fixedResiOffset=0,
            matchCount=-1
        ):
        rot, tran = self.getRotTrans(
            moving, 
            fixed,
            movingResiOffset=movingResiOffset,
            fixedResiOffset=fixedResiOffset,
            matchCount=matchCount
        )
        moving.transform(rot, tran)

    def getRotTrans(
            self, 
            moving, 
            fixed, 
            movingResiOffset=0, 
            fixedResiOffset=0,
            matchCount=-1
        ):
        movingChain = ElfinUtils.getChain(moving)
        ma = [
                al[0] for al in [[a for a in r.child_list if a.name == 'CA'] 
                for r in movingChain.child_list[movingResiOffset:(movingResiOffset+matchCount)]]
            ]

        fixedChain = ElfinUtils.getChain(fixed)
        fa = [
                al[0] for al in [[a for a in r.child_list if a.name == 'CA'] 
                for r in fixedChain.child_list[fixedResiOffset:(fixedResiOffset+matchCount)]]
            ]

        self.si.set_atoms(fa, ma)

        #   The rotation from BioPython is the
        # second dot operand instead of the 
        # conventional first dot operand.
        #   This means instead of R*v + T, the actual
        # transform is done with v'*R + T
        #   This is important to understand why I did
        # the rotation maths this way in the C++ GA
        return self.si.rotran

    def getRadii(self, pose):
        # Warning: this function assumes pose is centered!

        natoms = 0;
        rgSum = 0;
        maxCA = 0;

        nHeavy = 0;
        maxHeavy = 0;
        for a in pose.get_atoms():
            dist = np.linalg.norm(
                a.get_coord().astype('float64'));

            rgSum += dist;

            if(a.name =='CA'):
                maxCA = max(maxCA, dist);

            if(a.element != 'H'):
                maxHeavy = max(maxHeavy, dist);
                nHeavy = nHeavy + 1;

            natoms = natoms + 1;

        rg = rgSum / natoms;
        return OrderedDict([
            ('avgAll', rg),
            ('maxCA', maxCA),
            ('maxHeavy', maxHeavy)
        ]);

    def processSingle(self, filename):
        singleName = filename.split('/')[-1].replace('_0001.pdb', '')
        single = ElfinUtils.readPdb(filename)
        self.moveToOrigin(single)
        ElfinUtils.savePdb(
            single, 
            self.outputAlignedDir + '/single/' + singleName + '.pdb'
        )
        self.singlePDBs[singleName] = single

    def processPair(self, filename):
        # Step 1: Load pair and single structures
        pairName = filename.split('/')[-1].split('.')[0] \
            .replace('_0001', '')
        pair = ElfinUtils.readPdb(filename)

        singleNameA, singleNameB = pairName.split('-')
        singleA = self.singlePDBs[singleNameA]
        singleB = self.singlePDBs[singleNameB]

        rcA = ElfinUtils.getResidueCount(singleA)
        rcB = ElfinUtils.getResidueCount(singleB)
        rcPair = ElfinUtils.getResidueCount(pair)

        startResi = ElfinUtils.intFloor(float(rcA)/2)
        rcBEnd = ElfinUtils.intCeil(float(rcB)/2)
        endResi = rcPair - rcBEnd

        #   The fusionCount is the number of residues we use
        # to align pair to singleA. The high this number is, 
        # the more global our alignment is, which causes bad
        # disconnections in the chain. This is because in 
        # Synth we're inevitably fusing atom positions from 
        # different pairs into the same chain. Different pairs
        # have their single components stuck together using
        # and interface, the participation of which causes 
        # atom positions in the single components to differ 
        # from that of the original single module. 
        #   When we fuse different pairs together, each pair
        # is cut at 25% and 75% of their sequence in order to
        # be as far way to interfaces (0%, 50%, 100%) as 
        # possible. 
        #   The fusion alignment here is about aligning sub-
        # sequent pairs using a few residues before the 25% 
        # mark. The lower the fusionCount is, the fewer resi-
        # dues we use to align and the more local the align-
        # ment. However, if this number is too low the align-
        # ment could cause subsequent modules to overlap. 
        #   Through some experients I found that using 1/8 of
        # the length of singleB is a good balance between not
        # causing discontinuities and also not creating atom
        # overlaps.
        fusionCount = ElfinUtils.intCeil(float(rcB)/8)

        # Step 2: Move pair to align with first single
        #   This aligns pair by superimposing pair[0] 
        # with singleA
        #   We only want to align to the second quardrant 
        # of singleA's atoms. This is to be consistent
        # with step 5
        self.align(
            pair, 
            singleA, 
            matchCount=startResi
        )

        # Step 3: Get COM of the singleB as seen in the pair
        #    Only align the second quardrant of singleB
        # in order to be consistent with step 5
        comB = self.getCOM(
            singleB, 
            mother=pair, 
            childResiOffset=rcBEnd - fusionCount,
            motherResiOffset=rcA + rcBEnd - fusionCount,
            matchCount=fusionCount
        )

        # Step 4: Get radius for collision checks later:
        #           1. Avg dist to com (gyradius aka RG)
        #           2. Max dist from CA to com
        #           3. Max dist from any heavy stom (not H) to COM
        radA = self.getRadii(singleA)

        # Step 5: Get transformation of pair to the second single
        #   Pair is already aligned to first single so
        # there is no need for the first transformation
        #   You can check this is true by varifying that
        # self.getRotTrans(pair, singleA) has identity 
        # rotation and zero translation.
        #   Only align the second quardrant of singleB
        # in order to be consistent with the Synth 
        # script, where pairs are fused together by 
        # chopping the first and last quardrant of a 
        # pair. This means the second half of singleB 
        # is chopped off during fusion, while the 
        # first quardrant of singleB participates in 
        # interfacing. Therefore we align by uperimposing 
        # just the second quardrant
        rot, tran = self.getRotTrans(
            pair, 
            singleB, 
            movingResiOffset=(rcA + rcBEnd - fusionCount),
            fixedResiOffset=rcBEnd - fusionCount,
            matchCount=fusionCount
        )

        # Step 6: Save the aligned molecules
        #   Here the PDB format adds some slight 
        # floating point error. It is really old
        # and we should consider using mmCIF
        ElfinUtils.savePdb(
            pair, 
            self.outputAlignedDir + '/pair/' + pairName + '.pdb'
        )

        data = OrderedDict([
            ('comB',  comB.tolist()),
            ('rot',   rot.tolist()),
            ('tran',  tran.tolist())
        ])

        entry = self.pairsData.get(singleNameA, None)
        if entry == None:
            self.pairsData[singleNameA] = {}
            entry = self.pairsData.get(singleNameA, None)

        entry[singleNameB] = data

        singleDataA = self.singlesData.get(
            singleNameA,
             OrderedDict([
                ('linkCount', 0),
                ('radii', radA)
            ])
        );
        singleDataA['linkCount'] = singleDataA['linkCount'] + 1;
        self.singlesData[singleNameA] = singleDataA;

    def dumpXDB(self):
        toDump = OrderedDict([
            ('singlesData',  self.singlesData),
            ('pairsData',   self.pairsData)
            ])

        json.dump(toDump,
            open(self.outFile, 'w'),
            separators=(',', ':'),
            ensure_ascii=False,
            indent=4)

    def genFakeHubs(self):


    def dumpFakeHubs(self):
        toDump = OrderedDict([
            ('singlesData',  self.fakeHubsData.singles),
            ('pairsData',   self.fakeHubsData.pairs)
            ])

        json.dump(toDump,
            open(self.outFile, 'w'),
            separators=(',', ':'),
            ensure_ascii=False,
            indent=4)

    def run(self):
        if self.genFakeHubs:
            print '[MSG] Generating fake hubs'
            self.genFakeHubs()
            self.dumpFakeHubs()
        else:
            # Center all single modules
            files = glob.glob(self.singleDir + '/*.pdb')
            nSingles = len(files)

            for i in range(0, nSingles):
                print '[MSG] Centering single #{}/{}: {}' \
                    .format(i+1, nSingles, files[i])
                self.processSingle(files[i])

            # _0001 stands for relaxed PDBs
            files = glob.glob(self.pairDir + '/*.pdb')
            nPairs = len(files)
            for i in range(0, nPairs):
                print '[MSG] Aligning pair #{}/{}: {}' \
                    .format(i+1, nPairs, files[i])
                self.processPair(files[i])

            print '[Summary] {} singles, {} pairs'.format(nSingles, nPairs)

            self.dumpXDB()

if __name__ =='__main__': main()