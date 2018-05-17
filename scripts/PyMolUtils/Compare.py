
import numpy as np
from pymol import cmd

# Not sure how to just figure out where elfin is located
# So we need to load our library this way
elfinDir = '/Users/joy/src/elfin/'
elfinPyLibDir = elfinDir + '/src/python/'
elfinMovieDir = elfinDir + '/movieOutput/'
import imp
utils = imp.load_source('utils', elfinPyLibDir + '/utils.py')
Kabsch = imp.load_source('Kabsch', elfinPyLibDir + '/Kabsch.py')

def compare_sol(specFile, solCsv):
    if specFile.rfind('.csv') != -1:
        specPts = utils.readCsvPoints(specFile)
    elif specFile.rfind('.json') != -1:
        with open(specFile, 'r') as file:
            specPts = np.asarray(json.load(file)['coms'])
    else:
        print 'Unknown spec file format'

    solPts = utils.readCsvPoints(solCsv)

    # Centre both pts
    centredSpec = specPts - np.mean(specPts, axis=0)
    centredSol = solPts - np.mean(solPts, axis=0)

    # Draw specification
    draw_pts(centredSpec, color=[0.7,0,0])

    # Equalise sample points
    specUpPts = utils.upsample(centredSpec, centredSol)

    draw_pts(specUpPts, color=[0.5,0.5,0])

    # Find Kabsch rotation for solution -> spec
    R = Kabsch.kabsch(centredSpec, specUpPts)

    centredSpecR = np.dot(centredSpec, R)

    draw_pts(centredSpecR, color=[0,0.5,0.7])

    cmd.reset()
    cmd.set("depth_cue", 0)


cmd.extend("compare_sol", compare_sol)