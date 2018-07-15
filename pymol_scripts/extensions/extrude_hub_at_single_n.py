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
    def extrude_hub_at_single_n(single_name=None, hub_name=None, component_id=None):
        '''
        Extrudes a hub at the n-terminus of a single module.

        Args:
        - single_name - string
        - hub_name - string
        - component_id - string id of the module component inside the hub to
            extend into
        '''
        if single_name is None or \
            hub_name is None or \
            component_id is None:
            print(extrude_hub_at_single_n.__doc__)
        else:
            double_name = '-'.join([single_name, single_name])

            pdb_dir = os.getcwd() + '/../../resources/pdb_aligned/'

            cmd.load(pdb_dir + '/singles/' + single_name + '.pdb')
            cmd.set_name(single_name, 'single')
            cmd.load(pdb_dir + '/doubles/' + double_name + '.pdb')
            cmd.set_name(double_name, 'double')
            cmd.load(pdb_dir + '/doubles/' + double_name + '.pdb')
            cmd.set_name(double_name, 'double-o')
            cmd.load(pdb_dir + '/hubs/' + hub_name + '.pdb')
            cmd.set_name(hub_name, 'hub')

            xdb=utilities.read_json(os.getcwd() + '/../../resources/xdb.json')
            double_info = xdb['double_data'][single_name][single_name]

            # first, drop the double for reference
            tx('double', rot=double_info['rot'], tran_after=double_info['tran'])

            hub_comp_info = xdb['hub_data'][hub_name]['component_info']
            comp_a_cc = hub_comp_info[component_id]['c_connections'][single_name]

            # The frame drop order is important here

            # drop hub into component frame
            tx('hub', rot=comp_a_cc['rot'], tran_after=comp_a_cc['tran'])

            # drop hub into double's frame
            tx('hub', rot=double_info['rot'], tran_after=double_info['tran'])

            noclip()

            print('Extruded Hub {} Component {} at Single {}\'s N-Term'.\
                format(hub_name, component_id, single_name))

    @cmd.extend
    def extrude_hub_at_single_n_example():
        extrude_hub_at_single_n(single_name='D79', hub_name='D79_aC2_04', component_id='B')

    print('Extrude Hub At Single N Loaded')