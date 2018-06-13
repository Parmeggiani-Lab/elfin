# A PyMol mass export script
from pymol import cmd

import glob

def elfin_mass_obj_convert(srcDir, dstDir, ext='obj'):
    files = glob.glob(srcDir + '/single/*.pdb') + glob.glob(srcDir + '/pair/*.pdb') 

    cmd.set('auto_show_nonbonded', 'off')
    cmd.set('auto_show_selections', 'off')  
    cmd.set('auto_show_spheres', 'off')   
    cmd.set('auto_show_classified', 'off')  
    cmd.set('auto_show_lines', 'off')     

    fnInfo = []
    for f in files:
        cmd.load(f)
        name = '.'.join((re.split(r'/|\\', f)[-1]).split('.')[:-1])
        fnInfo.append((f, name))

    cmd.disable('all')
    cmd.show('cartoon')
    for fn in fnInfo:
        cmd.disable('all')
        cmd.enable(fn[1])
        cmd.save(fn[0].replace(srcDir, dstDir).replace('.pdb', '.' + ext))
        
    # Clear workspace
    cmd.delete('all')

# Clear workspace
cmd.reinitialize()

cmd.extend("elfin_mass_obj_convert", elfin_mass_obj_convert)