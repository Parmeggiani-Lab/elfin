#!/usr/bin/env python

import Bio.PDB
import ElfinUtils
import glob
import numpy as np
import codecs, json
import argparse
from collections import OrderedDict

def main():
    ap = argparse.ArgumentParser(description='Generates the xDB database from preprocessed single and double modules.');
    ap.add_argument('--outputAlignedDir', default='./resources/pdb_aligned/')
    ap.add_argument('--doublesDir', default='./resources/pdb_relaxed/doubles/')
    ap.add_argument('--singlesDir', default='./resources/pdb_relaxed/singles/')
    ap.add_argument('--output', default='./resources/xDB.json')
    args = ap.parse_args()

    xdbg = XDBGenrator(
        args.doublesDir, 
        args.singlesDir, 
        args.outputAlignedDir, 
        args.output)
    ElfinUtils.safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                doublesDir,
                singlesDir,
                outputAlignedDir,
                outFile):
        self.doublesDir     = doublesDir
        self.singlesDir     = singlesDir
        ElfinUtils.mkdir(outputAlignedDir)
        ElfinUtils.mkdir(outputAlignedDir     + '/doubles/')
        ElfinUtils.mkdir(outputAlignedDir     + '/singles/')
        self.outputAlignedDir  = outputAlignedDir
        self.outFile        = outFile
        self.si             = Bio.PDB.Superimposer()
        self.doublesData      = {}
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
            # This is for finding COM of a single inside a double
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

        #   The rotation from BioPython is the second dot operand instead of
        # the conventional first dot operand.
        #
        #   This means instead of R*v + T, the actual transform is done with
        # v'*R + T
        #
        #   This is important to understand why I did the rotation maths this
        # way in the C++ GA
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
        singleName = filename.split('/')[-1].replace('.pdb', '')
        single = ElfinUtils.readPdb(filename)
        self.moveToOrigin(single)
        ElfinUtils.savePdb(
            single, 
            self.outputAlignedDir + '/singles/' + singleName + '.pdb'
        )
        self.singlePDBs[singleName] = single

    def processDouble(self, filename):
        # Step 1: Load structures
        doubleName = filename.split('/')[-1].replace('.pdb', '') # _0001 should have been replaced befor thie step by some Shell script
        double = ElfinUtils.readPdb(filename)

        singleNameA, singleNameB = doubleName.split('-')
        singleA = self.singlePDBs[singleNameA]
        singleB = self.singlePDBs[singleNameB]

        rcA = ElfinUtils.getResidueCount(singleA)
        rcB = ElfinUtils.getResidueCount(singleB)
        rcDouble = ElfinUtils.getResidueCount(double)

        startResi = ElfinUtils.intFloor(float(rcA)/2)
        rcBEnd = ElfinUtils.intCeil(float(rcB)/2)
        endResi = rcDouble - rcBEnd

        #   The fusionCount is the number of residues we use to align double
        # to singleA. The high this number is, the more global our alignment
        # is, which causes bad disconnections in the chain. This is because in
        # Synth we're inevitably fusing atom positions from different doubles
        # into the same chain. Different doubles have their single components
        # stuck together using and interface, the participation of which
        # causes atom positions in the single components to differ from that
        # of the original single module.
        #
        #   When we fuse different doubles together, each double is cut at 25%
        # and 75% of their sequence in order to be as far way to interfaces
        # (0%, 50%, 100%) as possible.
        #
        #   The fusion alignment here is about aligning sub- sequent doubles
        # using a few residues before the 25% mark. The lower the fusionCount
        # is, the fewer resi- dues we use to align and the more local the
        # align- ment. However, if this number is too low the align- ment
        # could cause subsequent modules to overlap.
        #
        #   Through some experients I found that using 1/8 of the length of
        # singleB is a good balance between not causing discontinuities and
        # also not creating atom overlaps.
        fusionCount = ElfinUtils.intCeil(float(rcB)/8)

        # Step 2: Move double to align with first single
        #
        #   This aligns double by superimposing double[0] with singleA
        #
        #   We only want to align to the second quardrant of singleA's atoms.
        # This is to be consistent with step 5
        self.align(
            double, 
            singleA, 
            matchCount=startResi
        )

        # Step 3: Get COM of the singleB as seen in the double
        #
        #    Only align the second quardrant of singleB in order to be
        # consistent with step 5
        comB = self.getCOM(
            singleB, 
            mother=double, 
            childResiOffset=rcBEnd - fusionCount,
            motherResiOffset=rcA + rcBEnd - fusionCount,
            matchCount=fusionCount
        )

        # Step 4: Get radius for collision checks later:
        #           1. Avg dist to com (gyradius aka RG)
        #           2. Max dist from CA to com
        #           3. Max dist from any heavy stom (not H) to COM
        radA = self.getRadii(singleA)

        # Step 5: Get transformation of double to the second single
        #
        #   Double is already aligned to first single so there is no need for
        # the first transformation
        #
        #   You can check this is true by varifying that
        # self.getRotTrans(double, singleA) has identity rotation and zero
        # translation.
        #
        #   Only align the second quardrant of singleB in order to be
        # consistent with the Synth script, where doubles are fused together
        # by chopping the first and last quardrant of a double. This means the
        # second half of singleB is chopped off during fusion, while the first
        # quardrant of singleB participates in interfacing. Therefore we align
        # by uperimposing just the second quardrant
        rot, tran = self.getRotTrans(
            double, 
            singleB, 
            movingResiOffset=(rcA + rcBEnd - fusionCount),
            fixedResiOffset=rcBEnd - fusionCount,
            matchCount=fusionCount
        )

        # Step 6: Save the aligned molecules
        #
        #   Here the PDB format adds some slight floating point error. It is
        # really old and we should consider using mmCIF
        ElfinUtils.savePdb(
            double, 
            self.outputAlignedDir + '/doubles/' + doubleName + '.pdb'
        )

        data = OrderedDict([
            ('comB',  comB.tolist()),
            ('rot',   rot.tolist()),
            ('tran',  tran.tolist())
        ])

        entry = self.doublesData.get(singleNameA, None)
        if entry == None:
            self.doublesData[singleNameA] = {}
            entry = self.doublesData.get(singleNameA, None)

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
            ('singlesData', self.singlesData),
            ('doublesData', self.doublesData)
            ])

        json.dump(toDump,
            open(self.outFile, 'w'),
            separators=(',', ':'),
            ensure_ascii=False,
            indent=4)

    def run(self):
        # Center all single modules
        singleFiles = glob.glob(self.singlesDir + '/*.pdb')
        nSingles = len(singleFiles)
        for i in range(0, nSingles):
            print 'Centering single [{}/{}] {}' \
                .format(i+1, nSingles, singleFiles[i])
            self.processSingle(singleFiles[i])

        # _0001 stands for relaxed PDBs
        doubleFiles = glob.glob(self.doublesDir + '/*.pdb')
        nDoubles = len(doubleFiles)
        for i in range(0, nDoubles):
            print 'Aligning double [{}/{}] {}' \
                .format(i+1, nDoubles, doubleFiles[i])
            self.processDouble(doubleFiles[i])

        print 'Total: {} singles, {} doubles'.format(nSingles, nDoubles)

        self.dumpXDB()

if __name__ =='__main__': main()