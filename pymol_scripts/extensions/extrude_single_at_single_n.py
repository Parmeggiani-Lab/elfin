#!/usr/bin/env python3

#
# A PyMol extension script 
#

def main():
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
    from pelfin import utilities
    import numpy as np
    import os

    @cmd.extend
    def extrude_single_at_single_n(single_name=None, ext_single_name=None):
        '''
        Extrudes a single at the n-terminus of a single module.

        Args:
        - single_name - string name of the fixed single
        - ext_single_name - string name of the extension single
        '''
        if single_name is None or \
            ext_single_name is None:
            print(extrude_single_at_single_n.__doc__)
        else:
            double_name = '-'.join([ext_single_name, single_name])

            pdb_dir = os.getcwd() + '/../../resources/pdb_aligned/'

            cmd.load(pdb_dir + '/doubles/' + double_name + '.pdb')
            cmd.set_name(double_name, 'double')
            cmd.load(pdb_dir + '/singles/' + single_name + '.pdb')
            cmd.set_name(single_name, 'single')
            cmd.load(pdb_dir + '/singles/' + ext_single_name + '.pdb')
            cmd.set_name(ext_single_name, 'single-ext')
            
            xdb=utilities.read_json(os.getcwd() + '/../../resources/xdb.json')
            double_info = xdb['double_data'][ext_single_name][single_name]

            # first, drop the double (into its A frame) for reference
            tx('double', rot=double_info['rot'], tran_after=double_info['tran'])

            # extrude N term - drop into double's A frame
            tx('single-ext', rot=double_info['rot'], tran_after=double_info['tran'])

            cmd.disable('single-*')
            cmd.enable('single-ext')

            noclip()

            print('Extruded Single {} at Single {}\'s N-Term'.\
            format(ext_single_name, single_name))

    @cmd.extend
    def extrude_single_at_single_n_example():
        extrude_single_at_single_n(single_name='D79_j1_D54', ext_single_name='D79')

    print('Extrude Single At Single N Loaded')
