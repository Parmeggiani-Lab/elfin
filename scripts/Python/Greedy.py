#!/usr/bin/env python

import time
import re
import numpy as np
import glob, sys

haveCmd = False
try:
    from pymol import cmd
    haveCmd = True

    def safely(func, *args):
        try:
            func(*args)
        except Exception as e:
            print 'WTF... {}'.format(e)
except ImportError:
    print 'Could not import pymol cmd. Not running as pymol plugin...'

# Not sure how to just figure out where elfin is located
# So we need to load our library this way
elfinDir = '/Users/joy/src/elfin/'
elfinPyLibDir = elfinDir + '/src/python/'
elfinMovieDir = elfinDir + '/movieOutput/'
import imp
ElfinUtils = imp.load_source('ElfinUtils', elfinPyLibDir + '/ElfinUtils.py')
Kabsch = imp.load_source('Kabsch', elfinPyLibDir + '/Kabsch.py')

def main(pymolArgs=None):
    if haveCmd:
        sys.argv = ['pymol_script'] + pymolArgs

    if(len(sys.argv) < 2):
        print './Greedy.py <specFile.{json|pdb}> s=scale l=len'
        exit()

    specFile = sys.argv[1]
    fileExt = specFile[specFile.rfind('.'):]
    inputName = specFile[(specFile.rfind('/')+1):specFile.rfind('.')]

    if haveCmd:
        ElfinUtils.mkdir(elfinMovieDir + '/' + inputName + '/')

    xDB = ElfinUtils.readJSON(elfinDir + 'res/xDB.json')
    if fileExt == '.json':
        spec = ElfinUtils.readJSON(specFile)
        targetLen = len(spec['nodes'])
    elif fileExt == '.csv':
        with open(specFile, 'r') as file:
            pts = [[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n')]

        spec = {'coms': pts}
        npts = np.asarray(pts)

        (avgD, minD, maxD) = ElfinUtils.getXDBStat(xDB)
        # Use total length/avgD as heuristic. avgD is average xDB pair distance
        targetLen = int(np.ceil(sum(np.linalg.norm((npts-np.roll(npts, 1, axis=0))[1:], axis=1)) / avgD))
    else:
        print 'Unknown spec file type: {}'.format(fileExt)
        exit()

    # default options
    scale = 1.0
    userLen = -1
    for i in xrange(2, len(sys.argv)):
        arg = sys.argv[i]
        if arg.startswith('s='):
            scale = float(arg[(arg.rfind('=')+1):])
        elif arg.startswith('l='):
            userLen = int(arg[(arg.rfind('=')+1):])
        else:
            print 'Unrecognised flag: \"{}\"'.format(arg)
            exit()

    ElfinUtils.die(scale < 0, 'Scale must be > 0.0')

    # User specified length should overwrite 
    if userLen != -1:
        ElfinUtils.die(userLen < 3, 'Length must be > 3')
        targetLen = userLen
    else:
        targetLen = int(round(targetLen * scale))

    pairsDir    = elfinDir + 'res/centered_pdb/pair'
    singlesDir  = elfinDir + 'res/centered_pdb/single'
    designer    = GreedyDesigner(xDB, 'maxHeavy')
    startTime   = time.clock()
    (nodes,shape,score, fRot) = \
        designer.design(inputName,
                    spec, 
                    targetLen, 
                    scaleFactor=scale, 
                    singlesDir=singlesDir)
    print '{:.2f}s, score: {}'.format(time.clock() - startTime, score)

    if not haveCmd:
        ElfinUtils.makePdbFromNodes(xDB, nodes, pairsDir,
            specFile.replace(fileExt, ElfinUtils.suffixPdb(
                designer.__class__.__name__, 
                'Main',
                scale, 
                targetLen)),
            fRot)
            
class GreedyDesigner():
    frameCount = 0

    def __init__(self, xDB, collisionMeasure):
    	self.xDB                = xDB
        self.collisionMeasure   = collisionMeasure
        assert collisionMeasure in ElfinUtils.RadiiTypes

    def grow(self, lastNode, shape, newNode):
        rel = self.xDB['pairsData'][lastNode][newNode]
        shape = np.dot(shape, np.asarray(rel['rot'])) + rel['tran']
        shape = np.append(shape, [[0,0,0]], axis=0)

        return shape

    def transformPymolShape(self, idList, preTran, rot, postTran):
        ttt = np.concatenate(
                (np.concatenate((np.transpose(rot), 
                            np.transpose([
                                np.asarray(postTran).ravel().tolist()
                                ])), axis=1),
                [np.asarray(preTran).ravel().tolist() + [1]]));

        for single in idList:
            safely(
                cmd.transform_selection,
                single, ttt.ravel().tolist())

    def growShape(self, lastNode, idList, newNode):
        rel = self.xDB['pairsData'][lastNode][newNode]
        self.transformPymolShape(idList, [0,0,0], rel['rot'], rel['tran'])

        copyNodeName = newNode + '#' + str(len(idList))
        cmd.copy(copyNodeName, newNode)

        idList = np.append(idList, copyNodeName)

        return idList, rel['rot'], rel['tran'], copyNodeName

    def upsampleTarget(self, target):
        # Distances avg: 38.0111371052, min: 17.6585177454, max: 50.4035189874
        avgD = 38.0111371052
        minD = 17.6585177454
        maxD = 50.4035189874

        result = []
        for i in xrange(1, len(target)):
            result.append(target[i-1])

            delta = target[i]-target[i-1]
            dist = np.linalg.norm(delta)

            if dist < maxD:
                continue

            # Use min or avg!?
            n = int(np.floor(dist / avgD))

            partDelta = delta / n
            for j in xrange(1, n):
                result.append(target[i-1] + (j * partDelta))

        result.append(target[-1])
        return np.asarray(result)

    def evalShape(self, candidate, target, allowPerp=True):
        sumScore = 0
        for point in candidate:
            sumScore += ElfinUtils.minDistFromLine(point, target, allowPerp=allowPerp)
        return sumScore

    def loadSinglePdbsIntoPymol(self, singlesDir):
        # Load all singles
        ElfinUtils.die(singlesDir is None, 
            'Must provide singles directory when attempting to display in PyMol!')
        for singlePdb in glob.glob(singlesDir + '/*.pdb'):
            cmd.load(singlePdb)
        cmd.hide('everything', 'all')
        cmd.reset()

    def alignPymolView(self):
        view = (0.7483062148094177, 0.2898084819316864, 0.596701443195343, -0.0051952870562672615, 0.9020542502403259, -0.4315955936908722, -0.6633365750312805, 0.31986361742019653, 0.6765128374099731, 2.121180295944214e-05, 1.3768672943115234e-05, -1341.2861328125, 5.958309650421143, 15.377988815307617, 5.510014057159424, 1002.2947387695312, 1680.2774658203125, -20.0)
        cmd.set_view(view)

    def snapshotPymolShape(self, 
                        inputName,
                        idList, 
                        rot, 
                        tran, 
                        extension=None):

        if extension is not None:
            # Temporarily extend shape
            lastNodeName = idList[-1]
            lastNodeName = lastNodeName[:lastNodeName.rfind('#')]
            idList,extRot,extTran,extName = \
                self.growShape(lastNodeName, 
                    idList, 
                    extension)

        self.transformPymolShape(idList, tran, rot, [0, 0, 0])
        for single in idList:
            cmd.show('cartoon', single)
        self.alignPymolView()
        cmd.png(elfinMovieDir + inputName + '/frame' + str(self.frameCount))
        self.frameCount += 1

        # Restore shape so next transform will be correct
        self.transformPymolShape(idList, 
            [0,0,0], 
            np.linalg.inv(rot), 
            -1 * np.asarray(tran))
        for single in idList:
            cmd.hide('cartoon', single)

        if extension is not None:
            # Restore temporary extension
            idList = idList[:-1]
            cmd.delete(extName)
            self.transformPymolShape(idList, 
                                -1 * np.asarray(extTran), 
                                np.linalg.inv(extRot), 
                                [0,0,0])

    def design(self, 
            inputName, 
            spec,
            targetLen, 
            scaleFactor=1.0, 
            singlesDir=None):
        print 'Greedy design: target length {}, scale {}'.format(
            targetLen, scaleFactor)

        minTotalScore = float('inf')
        minScoreShape = None
        minScoreNodes = None

        # Load target and shift it by the first CoM 
        # so we can calculate score progressively
        target = spec.get('coms', None)
        assert target is not None
        target = np.asarray(target)
        target = target - target[0]
        target *= scaleFactor

        target = self.upsampleTarget(target)

        assert len(target) >= 2

        # Display in pymol if running as plugin
        if haveCmd:
            self.loadSinglePdbsIntoPymol(singlesDir);
            draw_axis()
            draw_pts(target)
            self.alignPymolView()

        # Try all starting pairs - we don't have information
        # to know which pair is best to start with
        pd = self.xDB['pairsData']
        candidate = np.asarray([[0, 0, 0]])

        bestStartPair = None
        bestStartScore = float('inf')
        bestStartKR = None
        for s1 in pd.keys():
            for s2 in pd[s1].keys(): 
                rel = pd[s1][s2]
                tmpShape = np.asarray([[0, 0, 0], rel['comB']])
                kR = Kabsch.kabsch(tmpShape, target[:2])
                tmpShape = np.dot(tmpShape, kR)
                s = self.evalShape(tmpShape, target, allowPerp=False)
                if s < bestStartScore:
                    bestStartPair = (s1, s2)
                    bestStartScore = s
                    bestStartKR = kR

        assert bestStartPair is not None

        if haveCmd:
            firstSingle = bestStartPair[0] + '#0'
            cmd.copy(firstSingle, bestStartPair[0])
            pymolCandidate,_,_,_ = self.growShape(bestStartPair[0],
                            [firstSingle], 
                            bestStartPair[1])
            self.snapshotPymolShape(inputName, 
                                pymolCandidate, 
                                bestStartKR, 
                                [0,0,0])

        candidate = self.grow(bestStartPair[0], [[0, 0, 0]], bestStartPair[1])
        nodes = [bestStartPair[0], bestStartPair[1]]
       

        # Main greedy construction loop
        totalScore = 0.0
        for i in xrange(2, targetLen):
            s1 = nodes[-1]

            bestNextNode = None
            bestNextScore = float('inf')
            atLeastOneNonColliding = False
            bestKR = None
            for s2 in pd[s1].keys(): 
                if ElfinUtils.checkCollision(self.xDB, 
                    self.collisionMeasure, 
                    nodes, 
                    s2, 
                    candidate):
                    continue

                atLeastOneNonColliding = True

                tmpShape = self.grow(s1, np.copy(candidate), s2)
                tmpShapeOrigin = tmpShape[0]
                tmpShape = tmpShape - tmpShapeOrigin

                matchLen = min(len(tmpShape), len(target))
                kR = Kabsch.kabsch(tmpShape[:matchLen], target[:matchLen])
                tmpShape = np.dot(tmpShape, kR)

                if haveCmd:
                    self.snapshotPymolShape(inputName,
                                        pymolCandidate, 
                                        kR, 
                                        -tmpShapeOrigin, 
                                        s2)

                s = self.evalShape(tmpShape, target, allowPerp=True)
                if s < bestNextScore:
                    bestNextNode = s2
                    bestNextScore = s
                    bestKR = kR

            if not atLeastOneNonColliding:
                print 'All next nodes collide! Cannot continue... Current length: {}'.format(i+1)
                print nodes
                print 'Available next nodes: {}'.format(pd[s2].keys())
                break

            assert bestNextNode is not None
            
            nodes.append(bestNextNode)
            candidate = self.grow(s1, candidate, bestNextNode)

            if haveCmd:
                pymolCandidate,_,_,_ = self.growShape(s1, 
                                                    pymolCandidate, 
                                                    bestNextNode)

            totalScore += bestNextScore

        # ElfinUtils.pauseCode()
        matchLen = min(len(candidate), len(target))
        kR = Kabsch.kabsch(candidate[:matchLen]-candidate[0], target[:matchLen])
        candidate = np.dot(candidate, kR)
        return nodes, candidate, totalScore, kR

if __name__ =='__main__': 
    ElfinUtils.safeExec(main)