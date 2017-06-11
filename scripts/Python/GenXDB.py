#!/usr/bin/env python

import Bio.PDB
import utils
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
    utils.safeExec(xdbg.run)

class XDBGenrator:

    def __init__(self,
                pairDir,
                singleDir,
                alignedLibDir,
                outFile,
                permissive=0):
        self.pairDir        = pairDir
        self.singleDir      = singleDir
        utils.mkdir(alignedLibDir)
        utils.mkdir(alignedLibDir     + '/pair/')
        utils.mkdir(alignedLibDir     + '/single/')
        self.alignedLibDir  = alignedLibDir
        self.outFile        = outFile
        self.si             = Bio.PDB.Superimposer()
        self.pairsData      = {}
        self.singlesData    = {}

    def getCOM(
            self, 
            child, 
            mother=None, 
            childAtomOffset=0, 
            motherAtomOffset=0,
            matchRange=-1
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
                movingAtomOffset=childAtomOffset, 
                fixedAtomOffset=motherAtomOffset,
                matchRange=matchRange
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
            movingAtomOffset=0, 
            fixedAtomOffset=0,
            matchRange=-1
        ):
        rot, tran = self.getRotTrans(
            moving, 
            fixed,
            movingAtomOffset=movingAtomOffset,
            fixedAtomOffset=fixedAtomOffset,
            matchRange=matchRange
        )
        moving.transform(rot, tran)

    def getRotTrans(
            self, 
            moving, 
            fixed, 
            movingAtomOffset=0, 
            fixedAtomOffset=0,
            matchRange=-1
        ):
        # First push the generators to desired locations
        maGen = moving.get_atoms()
        for i in xrange(0, movingAtomOffset):
            try:
                maGen.next()
            except StopIteration:
                die(True, 'Moving PDB Atom Offset too large')    

        faGen = fixed.get_atoms()
        for i in xrange(0, fixedAtomOffset):
            try:
                faGen.next()
            except StopIteration:
                die(True, 'Fixed PDB Atom Offset too large')    

        # Then fill in the arrays until either 
        # of the generators run out
        ma = []
        fa = []
        while True:
            try:
                maNext = maGen.next()
                faNext = faGen.next()
            except StopIteration:
                break

            if matchRange == -1 or len(ma) < matchRange:
                ma.append(maNext)
                fa.append(faNext)

        self.si.set_atoms(fa, ma)

        # Import note:
        # The rotation from BioPython seems be the
        # second dot operand instead of the 
        # conventional first dot operand!
        #
        # This means instead of R*v + T, the actual
        # transform is done with v'*R + T
        #
        # This has important to understand why I did
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

    def processPDB(self, filename):
        # Step 0: Load pair and single structures
        pairName = filename.split('/')[-1].split('.')[0] \
            .replace('_0001', '')
        pair = utils.readPdb(pairName, filename)

        singleNameA, singleNameB = pairName.split('-')
        singleA = utils.readPdb(
            singleNameA, 
            self.singleDir + singleNameA + '.pdb'
        )
        singleB = utils.readPdb(
            singleNameB, 
            self.singleDir + singleNameB + '.pdb'
        )

        atomCountA = utils.getAtomCount(singleA)
        atomCountB = utils.getAtomCount(singleB)

        # Step 1: Center the corresponding singles
        self.moveToOrigin(singleA)
        self.moveToOrigin(singleB)

        # Step 2: Move pair to align with first single
        # Note: this aligns pair by superimposing pair[0] 
        #       with singleA
        # Note: we only want to align the first half of 
        #       singleA's atoms, because the second half
        #       take part in interfacing with singleB, 
        #       and could be distorted differently in 
        #       different pairs
        self.align(pair, singleA, matchRange=utils.intFloor(atomCountA/2))

        # Step 3: Get COM of the singleB as seen in the pair
        # Note: only align the second half of singleB's atoms
        #       for the same reason as Step 2
        comB = self.getCOM(
            singleB, 
            pair, 
            childAtomOffset=utils.intFloor(atomCountB/2),
            motherAtomOffset=atomCountA + utils.intFloor(atomCountB/2)
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
        #       has identity rotation and zero translation. Also,
        #       comB should be at the origin.
        # Note: again, only match second half of singleB's atoms
        rot, tran = self.getRotTrans(
            pair, 
            singleB, 
            movingAtomOffset=atomCountA + utils.intFloor(atomCountB/2),
            fixedAtomOffset=utils.intFloor(atomCountB/2),
        )

        # Step 6: Save the centred and aligned molecules
        # Note: here the PDB format adds some slight 
        #       floating point error. It is really old
        #       and we should consider using mmCIF
        utils.savePdb(
            singleA, 
            self.alignedLibDir + '/single/' + singleNameA + '.pdb'
        )
        utils.savePdb(
            singleB, 
            self.alignedLibDir + '/single/' + singleNameB + '.pdb'
        )
        utils.savePdb(
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
        # _0001 stands for relaxed PDBs
        files = glob.glob(self.pairDir + '/*.pdb')
        nFiles = len(files)
        for i in range(0, nFiles):
            print '[XDBG] Processing file #{}/{}: {}' \
                .format(i+1, nFiles, files[i])
            self.processPDB(files[i])

        self.complexity = 1
        for s in self.singlesData:
            self.complexity = self.complexity * \
                self.singlesData.get(s)['linkCount']

        print '[XDBG] Complexity: {}'.format(self.complexity)

        self.dumpJSON()

if __name__ =='__main__': main()