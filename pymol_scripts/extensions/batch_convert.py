#!/usr/bin/env python3

#
# A PyMol extension script for batch converting objects (originally intended
# to convert into .obj models).
#

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()

in_pymol = False
try:
    import pymol
    in_pymol = True
except ImportError as ie:
    main()


if in_pymol: 
    from pymol import cmd

    import glob

    from elfinpy.utilities import make_dir

    @cmd.extend
    def batch_convert_modules(src_dir=None, dst_dir=None, ext='obj'):
        """Batch convert Elfin protein module PDBs.

        Args:
        - src_dir - input PDB directory (one that contains sub_dirs like singles,
            doubles, hubs)
        - dst_dir - output PDB directory
        - ext - file extension supported by PyMol
        """

        if src_dir is None or dst_dir is None:
            print(batch_convert.__doc__)
        else:
            # Clear workspace
            cmd.reinitialize()

            sub_dirs = ['singles', 'doubles', 'hubs'] # don't think we need cappings

            for sd in sub_dirs:
                make_dir(dst_dir + '/' + sd)

            files = [f for flist in [glob.glob(src_dir + '/{}/*.pdb'.format(sd)) for sd in sub_dirs] for f in flist]

            cmd.set('auto_show_nonbonded', 'off')
            cmd.set('auto_show_selections', 'off')  
            cmd.set('auto_show_spheres', 'off')   
            cmd.set('auto_show_classified', 'off')  
            cmd.set('auto_show_lines', 'off')     

            fn_info = []
            for f in files:
                cmd.load(f)
                name = '.'.join((re.split(r'/|\\', f)[-1]).split('.')[:-1])
                fn_info.append((f, name))

            cmd.disable('all')
            cmd.show('cartoon')
            for (file_path, module_name) in fn_info:
                cmd.enable(module_name)
                cmd.save(file_path.replace(src_dir, dst_dir).replace('.pdb', '.' + ext))
                cmd.disable(module_name)
                
            # Clear workspace
            cmd.delete('all')

    print('Batch Convert Loaded')