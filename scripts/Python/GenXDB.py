#!/usr/bin/env python

import Bio.PDB
import ElfinUtils
import glob
import numpy as np
import codecs, json
from collections import OrderedDict

def main():
    pairDir         = './res/relaxed/pair/'
    singleDir       = './res/relaxed/single/'
    alignedLibDir   = './res/aligned/'
    outFile         = './res/xDB.json'
    xdbg = XDBGenrator(pairDir, singleDir, alignedLibDir, outFile)
    ElfinUtils.safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                alignedLibDir,
                outFile,
                permissive=0):
        self.pairDir        = pairDir
        self.singleDir      = singleDir
        ElfinUtils.mkdir(alignedLibDir)
        ElfinUtils.mkdir(alignedLibDir     + '/pair/')
        ElfinUtils.mkdir(alignedLibDir     + '/single/')
        self.alignedLibDir  = alignedLibDir
        self.outFile        = outFile
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
        # First push the generators to desired locations
        mrGen = moving.get_residues()
        for i in xrange(0, movingResiOffset):
            try:
                mrGen.next()
            except StopIteration:
                die(True, 'Moving residue offset too large')    

        frGen = fixed.get_residues()
        for i in xrange(0, fixedResiOffset):
            try:
                frGen.next()
            except StopIteration:
                die(True, 'Fixed residue offset too large')    

        # Then fill in the arrays until either 
        # of the generators run out
        ma = []
        fa = []
        while True:
            if matchCount != -1 and len(ma) >= matchCount:
                break

            try:
                for mcas in [a for a in mrGen.next().get_atoms() if a.name == 'CA']:
                    ma.append(mcas)

                for fcas in [a for a in frGen.next().get_atoms() if a.name == 'CA']:
                    fa.append(fcas)
            except StopIteration:
                break

        self.si.set_atoms(fa, ma)

        # Import note:
        # The rotation from BioPython is the
        # second dot operand instead of the 
        # conventional first dot operand.
        #
        # This means instead of R*v + T, the actual
        # transform is done with v'*R + T
        #
        # This is important to understand why I did
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
        single = ElfinUtils.readPdb(singleName, filename)
        self.moveToOrigin(single)
        ElfinUtils.savePdb(
            single, 
            self.alignedLibDir + '/single/' + singleName + '.pdb'
        )
        self.singlePDBs[singleName] = single

    def processPair(self, filename):
        # Step 1: Load pair and single structures
        pairName = filename.split('/')[-1].split('.')[0] \
            .replace('_0001', '')
        pair = ElfinUtils.readPdb(pairName, filename)

        singleNameA, singleNameB = pairName.split('-')
        singleA = self.singlePDBs[singleNameA]
        singleB = self.singlePDBs[singleNameB]

        rcA = ElfinUtils.getResidueCount(singleA)
        rcB = ElfinUtils.getResidueCount(singleB)
        rcPair = ElfinUtils.getResidueCount(pair)

        startResi = ElfinUtils.intFloor(float(rcA)/2)
        rcBEnd = ElfinUtils.intCeil(float(rcB)/2)
        endResi = rcPair - rcBEnd

        # Note: fusionCount is not a magic number. 
        #   We used to align pairs to singleB by using the first 
        # half of singleB's residues. However, during Synth 
        # singleB's first half may also participate in an 
        # interface. 
        #   When singleB takes part in two interfaces (both c-
        # term and n-term) then we must use the most central res-
        # idues for alignment. The more resi-dues we use, the 
        # less tightly-fit the alignment becomes (since Kabsch 
        # tries to minimise global RMSD). Since 3 is the mininum 
        # number for 3D superposition, while still ensuring
        # that chains don't break, 3 is chosen here.
        #   Numbers like ElfinUtils.intFloor(float(rcB)/8) also 
        # seems to work, but there is not reason to believe under 
        # some designs that large a number would not introduce 
        # chain breakage
        fusionCount = 3

        # Step 2: Move pair to align with first single
        # Note: this aligns pair by superimposing pair[0] 
        #       with singleA
        # Note: we only want to align to the second quardrant 
        #       of singleA's atoms. This is to be consistent
        #       with step 5
        self.align(
            pair, 
            singleA, 
            matchCount=startResi
        )

        # Step 3: Get COM of the singleB as seen in the pair
        # Note: only align the second quardrant of singleB
        #       in order to be consistent with step 5
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
        # Note: pair is already aligned to first single so
        #       there is no need for the first transformation
        #       You can check this is true by varifying that
        #           self.getRotTrans(pair, singleA)
        #       has identity rotation and zero translation.
        # Note: only align the second quardrant of singleB
        #       in order to be consistent with the Synth 
        #       script, where pairs are fused together by 
        #       chopping the first and last quardrant of a 
        #       pair. This means the second half of singleB 
        #       is chopped off during fusion, while the 
        #       first quardrant of singleB participates in 
        #       interfacing. Therefore we align by
        #       superimposing just the second quardrant
        rot, tran = self.getRotTrans(
            pair, 
            singleB, 
            movingResiOffset=(rcA + rcBEnd - fusionCount),
            fixedResiOffset=rcBEnd - fusionCount,
            matchCount=fusionCount
        )

        # Step 6: Save the aligned molecules
        # Note: here the PDB format adds some slight 
        #       floating point error. It is really old
        #       and we should consider using mmCIF
        # Note: mmCIF still uses 3 decimals places; 
        #       perhaps we should to find a way to 
        #       increase that precision
        ElfinUtils.savePdb(
            pair, 
            self.alignedLibDir + '/pair/' + pairName + '.pdb'
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

        if(singleNameB != singleNameA):
            singleDataB = self.singlesData.get(singleNameB,
                 OrderedDict([
                    ('linkCount', 0),
                    ('radii', radA)
                    ]));
            self.singlesData[singleNameB] = singleDataB;

        # interact(globals(), locals())

    def dumpJSON(self):
        toDump = OrderedDict([
            ('singlesData',  self.singlesData),
            ('complexity',  self.complexity),
            ('pairsData',   self.pairsData)
            ])

        json.dump(toDump,
            open(self.outFile, 'w'),
            separators=(',', ':'),
            ensure_ascii=False,
            indent=4)

    def run(self):
        # Center all single modules
        files = glob.glob(self.singleDir + '/*.pdb')
        nFiles = len(files)
        for i in range(0, nFiles):
            print '[XDBG] Centering single #{}/{}: {}' \
                .format(i+1, nFiles, files[i])
            self.processSingle(files[i])

        # _0001 stands for relaxed PDBs
        files = glob.glob(self.pairDir + '/*.pdb')
        nFiles = len(files)
        for i in range(0, nFiles):
            print '[XDBG] Aligning pair #{}/{}: {}' \
                .format(i+1, nFiles, files[i])
            self.processPair(files[i])

        self.complexity = 1
        for s in self.singlesData:
            self.complexity = self.complexity * \
                self.singlesData.get(s)['linkCount']

        print '[XDBG] Complexity: {}'.format(self.complexity)

        self.dumpJSON()

if __name__ =='__main__': main()