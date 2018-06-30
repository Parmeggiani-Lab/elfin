
from pymol import cmd

import glob

from pyelfin import make_dir

def elfin_mass_obj_convert(src_dir, dst_dir, ext='obj'):
  subdirs = ['singles', 'doubles', 'hubs'] # don't think we need cappings

  for sd in subdirs:
    make_dir(dst_dir + '/' + sd)

  files = [f for flist in [glob.glob(src_dir + '/{}/*.pdb'.format(sd)) for sd in subdirs] for f in flist]

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
  for (file_path, module_name) in fnInfo:
    cmd.enable(module_name)
    cmd.save(file_path.replace(src_dir, dst_dir).replace('.pdb', '.' + ext))
    cmd.disable(module_name)
    
  # Clear workspace
  cmd.delete('all')

# Clear workspace
cmd.reinitialize()

cmd.extend("elfin_mass_obj_convert", elfin_mass_obj_convert)