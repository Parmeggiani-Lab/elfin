#!/usr/bin/env python3

#
# A PyMol extension script to test extrusion of a hub from a single module's
# c-term
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
    from elfinpy import utilities
    import os

    @cmd.extend
    def extrude_hub_at_single_c(single_name=None, hub_name=None, component_id=None):
        """Extrudes a hub at the c-terminus of a single module.

        Args:
        - single_name - string
        - hub_name - string
        - component_id - string, indicating which module component inside the hub
            to extend into
        """
        if single_name is None or \
            hub_name is None or \
            component_id is None:
            print(extrude_hub_at_single_c.__doc__)
        else:
            double_name = '-'.join([single_name, single_name])

            pdb_dir = os.getcwd() + '/../../resources/pdb_aligned/'

            cmd.load(pdb_dir + '/singles/' + single_name + '.pdb')
            cmd.set_name(single_name, 'single')
            cmd.load(pdb_dir + '/doubles/' + double_name + '.pdb')
            cmd.set_name(double_name, 'double')
            cmd.load(pdb_dir + '/hubs/' + hub_name + '.pdb')
            cmd.set_name(hub_name, 'hub')
            
            xdb=utilities.read_json(os.getcwd() + '/../../resources/xdb.json')
            hub_comp_info = xdb['hub_data'][hub_name]['component_info']
            comp_a_cc = hub_comp_info[component_id]['n_connections'][single_name]

            tx('hub', rot=comp_a_cc['rot'], tran_after=comp_a_cc['tran'])

            noclip()

            print('Extruded Hub {} Component {} at Single {}\'s C-Term'.\
                format(hub_name, component_id, single_name))

    @cmd.extend
    def extrude_hub_at_single_c_example():
        extrude_hub_at_single_c(single_name='D4', hub_name='D4_C3_02', component_id='B')

    print('Extrude Hub At Single Loaded')