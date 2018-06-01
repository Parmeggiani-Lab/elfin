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

class ElfinGraph():
    """Representation of a chain. It might be tree-like or even cyclic"""
    def __init__(self, name='', nodes=[]):
        self.name = name
        self.nodes = nodes

    def transform(self, rot, tran):
        for n in self.nodes:
            n.transform(rot, tran)

class ElfinNode():
    """Representation of a single module instance and stores info about connectivity"""
    def __init__(
        self, 
        id, 
        name, 
        trim=(False, False),
        cap=None,
        ctermNodeId=-1, 
        rot=[[1,0,0],[0,1,0],[0,0,1]], 
        tran=[0,0,0],
        extraAttachments=[]
    ):
        self.id = id
        self.name = name
        self.trim = trim
        self.cap = cap
        self.ctermNodeId = ctermNodeId
        self.rot = rot
        self.tran = tran

        # Default is to cap the end that is not trimmed
        if self.cap == None:
            self.cap = (not trim[0], not trim[1])

        # Error checking
        if self.id < 0:
            raise ValueError('Node ID should never be negative: id={}'.format(self.id))

        if len(self.trim) != 2:
            raise ValueError('ElfinNode trim vector length != 2: trim={}'.format(self.trim))

        if not self.trim[0] and not self.trim[1]:
            print ('Warning: ElfinNode trim vector both ends are NOT trimmed. '
            'This should only happen if the chain has one single node, which '
            'is not thought to be common.')

        for i in xrange(0,2):
            if self.trim[i] and self.cap[i]:
                raise ValueError('Cannot cap a trimmed end[{}]: name={}, id={}'
                    .format(i, self.name, self.id))

    def transform(self, rot, tran):
        self.rot = (np.dot(self.rot, rot)).tolist()
        self.tran = (np.dot(self.tran, rot) + tran).tolist()

def genPymolTxm(rot, tran):
    rotTp = np.transpose(rot)
    rotTpTran = np.append(rotTp, np.transpose([tran]), axis=1)
    pymolRotMat = np.append(rotTpTran, [[0,0,0,1]], axis=0)
    return '[' + ', '.join(map(str, pymolRotMat.ravel())) + ']'


def getChain(struct, chainName='A'):
    return struct.child_list[0].child_dict[chainName]

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
    upsampled = np.zeros([0, 3])

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
#     upsampled = np.zeros([0, 3])

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

def readCsvPoints(csvFile):
    pts = []
    
    with open(csvFile, 'r') as file:
        pts = np.asarray([[float(n) for n in re.split(', *| *', l.strip())] for l in file.read().split('\n') if len(l) > 0])
    
    return pts

def readCsv(filePath, delim=','):
    rows = []
    with open(filePath) as csvfile:
        sreader = csv.reader(csvfile, delimiter=delim)
        for r in sreader:
            rows.append(r)

    return rows

def saveCsv(npArray, saveFile, delimiter=' '):
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

def readPdb(
    inFile,
    pdbName=None,
    permissive=0
):
    if pdbName == None:
        pdbName = inFile.split('/')[-1].replace('.', '_')
    parser = Bio.PDB.PDBParser(permissive)
    structure = parser.get_structure(pdbName, inFile)
    return structure

def saveCif(struct, saveFile):
    io = Bio.PDB.mmcifio.MMCIFIO()
    io.set_structure(struct)
    io.save(saveFile + ('' if saveFile.endswith('.cif') else '.cif'))

def savePdb(struct, saveFile):
    io = Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(saveFile)
    io.save(saveFile + ('' if saveFile.endswith('.pdb') else '.pdb'))

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
