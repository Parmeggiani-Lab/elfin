
import numpy as np
from pymol import cmd

from pyelfin import *

def compare_solutions(spec_file, sol_csv_file):
    if spec_file.rfind('.csv') != -1:
        spec_pts = utils.read_csv_points(spec_file)
    elif spec_file.rfind('.json') != -1:
        with open(spec_file, 'r') as file:
            spec_pts = np.asarray(json.load(file)['coms'])
    else:
        print 'Unknown spec file format'

    sol_pts = utils.read_csv_points(sol_csv_file)

    # Centre both pts
    centred_spec = spec_pts - np.mean(spec_pts, axis=0)
    centred_sol = sol_pts - np.mean(sol_pts, axis=0)

    # Draw specification
    draw_pts(centred_spec, color=[0.7,0,0])

    # Equalise sample points
    specUpPts = utils.upsample(centred_spec, centred_sol)

    draw_pts(specUpPts, color=[0.5,0.5,0])

    # Find Kabsch rotation for solution -> spec
    R = kabsch.run_kabsch(centred_spec, specUpPts)

    centredSpecR = np.dot(centred_spec, R)

    draw_pts(centredSpecR, color=[0,0.5,0.7])

    cmd.reset()
    cmd.set("depth_cue", 0)


cmd.extend("compare_solutions", compare_solutions)