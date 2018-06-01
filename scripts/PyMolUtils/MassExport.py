# A PyMol mass export script
from pymol import cmd

import glob

def mass_export(srcDir, dstDir, ext='obj'):
    files = glob.glob(srcDir + '/*.pdb')
    for f in files:
        cmd.load(f)

        name = '.'.join((re.split(r'/|\\', f)[-1]).split('.')[:-1])

        print 'f=' + f
        print 'name=' + name
        cmd.save(dstDir + '/' + name + '.' + ext)
        
        # Clear workspace
        cmd.delete('all')

# Clear workspace
cmd.reinitialize()

cmd.extend("mass_export", mass_export)