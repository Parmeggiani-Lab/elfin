import inspect, sys, code, traceback
import json
import os
import Bio.PDB
import numpy as np
import csv
import re
import argparse

RadiiTypes = ['avgAll', 'maxCA', 'maxHeavy']
INF = float('inf')

def getResidueCount(pdb):
    return sum([len(c.child_list) for c in pdb.child_list[0].child_list])

def getAtomCount(pdb):
    i = 0
    for a in pdb.get_atoms():
        i += 1
    return i

def intCeil(f):
    return int(np.ceil(f))

def intFloor(f):
    return int(np.floor(f))

def upsample(spec, pts):
    N = len(spec)

    morePoints, fewerPoints = (np.copy(spec), np.copy(pts))

    # Compute longer shape total length
    mpTotalLength = 0.0
    for i in xrange(1, N):
        mpTotalLength += np.linalg.norm(morePoints[i] - morePoints[i - 1])

    if mpTotalLength == INF:
        print 'Something fishy... mpTotalLength is inf!'

    fpTotalLength = 0.0
    for i in xrange(1, len(fewerPoints)):
        fpTotalLength += np.linalg.norm(fewerPoints[i] - fewerPoints[i - 1])

    if mpTotalLength == INF:
        print 'Something fishy... fpTotalLength is inf!'

    # Upsample fewerPoints
    upsampled = np.empty([0, 3])

    # First and last points are the same
    upsampled = np.append(upsampled, [fewerPoints[0]], axis=0)

    mpProportion = 0.0
    fpProportion = 0.0
    mpi = 1
    for i in xrange(1, len(fewerPoints)):
        baseFpPoint = fewerPoints[i - 1]
        nextFpPoint = fewerPoints[i]
        baseFpProportion = fpProportion
        fpSegment = np.linalg.norm(nextFpPoint - baseFpPoint) / fpTotalLength
        vec = nextFpPoint - baseFpPoint

        fpProportion += fpSegment
        while (mpProportion <= fpProportion and mpi < N):
            mpSegment = \
                np.linalg.norm(morePoints[mpi] - morePoints[mpi - 1]) \
                / mpTotalLength

            if (mpProportion + mpSegment) > fpProportion:
                break
            mpProportion += mpSegment

            s = (mpProportion - baseFpProportion) / fpSegment
            upsampled = np.append(upsampled, [baseFpPoint + (vec * s)], axis=0)

            mpi += 1

    # Sometimes the last node is automatically added
    if len(upsampled) < N:
        upsampled = np.append(upsampled, [fewerPoints[-1]], axis=0)

    return upsampled

# def upsample(pts1, pts2):
#     # Upsample the shape with fewer points
#     pts1IsLonger = len(pts1) > len(pts2)
#     N = max(len(pts1), len(pts2))

#     # Use a proportion based algorithm because
#     # we want to assume both shapes are roughly
#     # the same length, but not exactly
#     if pts1IsLonger:
#         morePoints, fewerPoints = (np.copy(pts1), np.copy(pts2))
#     else:
#         morePoints, fewerPoints = (np.copy(pts2), np.copy(pts1))

#     # Compute longer shape total length
#     mpTotalLength = 0.0
#     for i in xrange(1, N):
#         mpTotalLength += np.linalg.norm(morePoints[i] - morePoints[i - 1])

#     if mpTotalLength == INF:
#         print 'Something fishy... mpTotalLength is inf!'

#     fpTotalLength = 0.0
#     for i in xrange(1, len(fewerPoints)):
#         fpTotalLength += np.linalg.norm(fewerPoints[i] - fewerPoints[i - 1])

#     if mpTotalLength == INF:
#         print 'Something fishy... fpTotalLength is inf!'

#     # Upsample fewerPoints
#     upsampled = np.empty([0, 3])

#     # First and last points are the same
#     upsampled = np.append(upsampled, [fewerPoints[0]], axis=0)

#     mpProportion = 0.0
#     fpProportion = 0.0
#     mpi = 1
#     for i in xrange(1, len(fewerPoints)):
#         baseFpPoint = fewerPoints[i - 1]
#         nextFpPoint = fewerPoints[i]
#         baseFpProportion = fpProportion
#         fpSegment = np.linalg.norm(nextFpPoint - baseFpPoint) / fpTotalLength
#         vec = nextFpPoint - baseFpPoint

#         fpProportion += fpSegment
#         while (mpProportion <= fpProportion and mpi < N):
#             mpSegment = \
#                 np.linalg.norm(morePoints[mpi] - morePoints[mpi - 1]) \
#                 / mpTotalLength

#             if (mpProportion + mpSegment) > fpProportion:
#                 break
#             mpProportion += mpSegment

#             s = (mpProportion - baseFpProportion) / fpSegment
#             upsampled = np.append(upsampled, [baseFpPoint + (vec * s)], axis=0)

#             mpi += 1

#     # Sometimes the last node is automatically added
#     if len(upsampled) < N:
#         upsampled = np.append(upsampled, [fewerPoints[-1]], axis=0)

#     if pts1IsLonger:
#         pts2 = upsampled 
#     else: 
#         pts1 = upsampled

#     return pts1, pts2

def readCSVPoints(csvFile):
    pts = []
    
    with open(csvFile, 'r') as file:
        pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
    
    return pts

def saveCSV(npArray, saveFile, delimiter=' '):
    with open(saveFile, 'wb') as csvFile:
        wt = csv.writer(csvFile, delimiter=delimiter)
        for row in npArray:
            wt.writerow(row)

def canConvertToFloat(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

# Credits to http://stackoverflow.com/questions/2597278/python-load-variables-in-a-dict-into-namespace
class Bunch(object):
  def __init__(self, adict):
    self.__dict__.update(adict)

def floatApproximates(a, b, error=1e-6):
    return abs(a-b) < error

def realPath(path):
    return os.path.realpath(path)

def normPath(path):
    return os.path.normpath(path)

def suffixPdb(className, fromFunction, scale, targetLen):
    return '_{}_{}_s{}_l{}.pdb'.format(
            className, 
            fromFunction,
            scale, 
            targetLen)

def minDistFromLine(point, linePoints, allowPerp=True):
    minDist = float('inf')
    for i in xrange(1, len(linePoints)):
        lineSeg = (linePoints[i-1], linePoints[i])

        # First determine whether point is outside line segment regime
        v = lineSeg[1] - lineSeg[0]
        w = point - lineSeg[0]

        if allowPerp:
            c1 = np.dot(w, v)
            if c1 <= 0: # before lineSeg[0]
                dist = np.linalg.norm(w)
            else:
                c2 = np.dot(v, v)
                if c2 <= c1: # after lineSeg[1]
                    dist = np.linalg.norm(point - lineSeg[1])
                else:
                    # If not outside, then calculate perpendicular distance
                    b = c1 / c2
                    pol = lineSeg[0] + b*v
                    dist = np.linalg.norm(point - pol)
        else:
            dist = min(np.linalg.norm(point - lineSeg[0]), np.linalg.norm(point - lineSeg[1]))

        if dist < minDist:
            minDist = dist

    return minDist

def die(condition, str):
    if condition:
        print str
        exit()

def checkCollision(xdb, collisionMeasure, nodes, newNode, shape):
    newCOM = xdb['pairsData'][nodes[-1]][newNode]['comB']

    # previous node PAIR (not just single node!) is inherently non-colliding
    for i in xrange(0, len(nodes) - 2):
        comDist = np.linalg.norm(shape[i] - newCOM)
        collisionDist = xdb['singlesData'][newNode]['radii'][collisionMeasure] + \
                            xdb['singlesData'][nodes[i]]['radii'][collisionMeasure]

        if comDist < collisionDist:
            return True

    return False

def makePdbFromNodes(xdb, nodes, pairsDir, singlesDir, saveFile=None, fRot=None, movieMode=False):
    # Consturct a protein using single modules and xDB relationships

    # Use the first pdb model as an empty host for all single pdb chains
    motherName = nodes[0]
    pdbFile = singlesDir + '/' + motherName + '.pdb'
    motherPdb = readPdb(motherName, pdbFile)
    motherModel = motherPdb.child_list[0]
    motherModel.detach_child('A')
    assert(len(motherModel.child_list) == 0)

    motherChain = Bio.PDB.Chain.Chain('A')
    motherModel.add(motherChain)

    moviePdbs = []
    residueUid = 1

    comShape = np.empty([1, 3])
    startingPoint = np.zeros(3)

    si = Bio.PDB.Superimposer()

    chainLenDigits = len(str(len(nodes)))
    for i in xrange(0, len(nodes) - 1):
        currNode = nodes[i]
        nextNode = nodes[i+1]
        pairName = currNode + '-' + nextNode
        rel = xdb['pairsData'][currNode][nextNode]
        
        # Append new point at origin
        comShape = np.append(comShape, [[0,0,0]], axis=0)
        
        if movieMode:
            singlePdb = readPdb(
                currNode, 
                singlesDir + '/' + currNode + '.pdb'
            )
            moviePdbs.append(singlePdb)
            for pdb in moviePdbs:
            	pdb.transform(np.asarray(rel['rot']), rel['tran'])

            # Don't forget to add the last node
            if i == len(nodes) - 2:
                singlePdb = readPdb(
                    nextNode, 
                    singlesDir + '/' + nextNode + '.pdb'
                )
                moviePdbs.append(singlePdb)
        else:
            singleA = readPdb(
                currNode,
                singlesDir + '/' + currNode + '.pdb'
            )
            singleB = readPdb(
                nextNode,
                singlesDir + '/' + nextNode + '.pdb'
            )
            pair = readPdb(
                pairName,
                pairsDir + '/' + pairName + '.pdb'
            )

            resiCountA = getResidueCount(singleA)
            resiCountB = getResidueCount(singleB)
            resiCountPair = getResidueCount(pair)

            if i == 0:
                # First pair: ignore leading trim
                startResi = 0
                endResi = resiCountPair - intCeil(float(resiCountB)/2)
            elif i == len(nodes) - 2:
                # Last pair: ignore trailing trim
                startResi = intFloor(float(resiCountA)/2)
                endResi = resiCountPair
            else:
                # Trim half of singleA's residues and extend
                # to half of singleB's residues
                startResi = intFloor(float(resiCountA)/2)
                endResi = resiCountPair - intCeil(float(resiCountB)/2)

            pairChain = pair.child_list[0].child_dict['A']
            pairMidKeys = [r.id for r in pairChain.child_list[startResi:endResi]]
            pairMidRes = [r for r in pairChain.child_list[startResi:endResi]]

            if i > 0:
                # Extract atoms from last pair end
                prevCAs = [a for r in motherChain.child_list[-startResi:] for a in r if a.name == 'CA']
                currCAs = [a for r in pairChain.child_list[:startResi] for a in r if a.name == 'CA']

                # print '{}\n{}'.format([r.get_resname() for r in motherChain.child_list[-startResi:]], \
                #     [r.get_resname() for r in pairChain.child_list[:startResi]])
                # pauseCode()

                # Move mother slightly to align with current pair CAs
                si.set_atoms(currCAs, prevCAs)
                rot, tran = si.rotran
                motherPdb.transform(rot, tran)

            [pairChain.detach_child(k) for k in pairMidKeys]

            for r in pairMidRes:
                r.id = (r.id[0], residueUid, r.id[2]) 
                motherChain.add(r)
                residueUid += 1
            
            motherPdb.transform(np.asarray(rel['rot']), rel['tran'])

        # It seems sometimes np.empty() gives weirdly large values..
        # print comShape
        comShape = np.dot(comShape, np.asarray(rel['rot'])) + rel['tran']

        startingPoint = np.dot(startingPoint, np.asarray(rel['rot'])) + rel['tran']
        print 'Pair #{}:   {}'.format(str(i+1).ljust(chainLenDigits), pairName.ljust(16))
        # print 'From file: {}'.format(pairsDir + '/' + pairName + '.pdb')
        # print 'Residues: {} - {}'.format(startResi, endResi)

       
    if fRot is not None:
        if movieMode:
            motherPdb.transform(np.eye(3), -startingPoint)
            motherPdb.transform(np.asarray(fRot), np.zeros(3))
        else:
            for pdb in moviePdbs:
                pdb.transform(np.eye(3), -startingPoint)
                pdb.transform(np.asarray(fRot), np.zeros(3))

    if not movieMode:
        if saveFile is not None:
            savePdb(motherPdb, saveFile)
        return motherPdb, comShape
    else:
        if saveFile is not None:
            pdbId = 0
            saveFileDotIndex = saveFile.rfind('.')
            for pdb in moviePdbs:
                savePartFile = saveFile[0:saveFileDotIndex] + \
                    'part' + str(pdbId) + \
                    saveFile[saveFileDotIndex:]
                savePdb(pdb, savePartFile)
                pdbId = pdbId + 1
        return moviePdbs, comShape

def getXDBStat(xDB):
    # xdb = readJSON('res/xDB.json')

    pd = xDB['pairsData']
    dists = []
    for s1 in pd.keys():
        for s2 in pd[s1].keys():
            dists.append(np.linalg.norm(pd[s1][s2]['comB']))

    return np.average(dists), min(dists), max(dists)

def readJSON(filename):
    with open(filename, 'r') as openFile:
        return json.load(openFile)

def readPdb(customName, inFile, permissive=0):
    parser = Bio.PDB.PDBParser(permissive)
    structure = parser.get_structure(customName, inFile)
    return structure

def savePdb(struct, saveFile):
    io = Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(saveFile)

def interact(globalVars=None, localsVars=None):
    print "Entering interactive mode"
    print
    if(localsVars == None):
        localsVars = locals()
    if(globalVars == None):
        globalVars = globals()
    code.interact(local=dict(globalVars, **localsVars))

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def pauseCode(frame=None):
    print '---------pauseCode()---------'
    if frame is None:
        # Use current frame (one above the exception wrapper)
        frame = inspect.currentframe().f_back
    
    ns = dict(frame.f_globals)
    ns.update(frame.f_locals)
    code.interact(local=ns)

def safeExec(func, *args):
    try:
        func(*args)
    except Exception as e:
        print '---------safeExec() caught exception---------'

        # Find last (failed) inner frame
        type, value, tb = sys.exc_info()
        last_frame = lambda tb=tb: last_frame(tb.tb_next) if tb.tb_next else tb
        frame = last_frame().tb_frame
        traceback.print_exc()
        pauseCode(frame)
