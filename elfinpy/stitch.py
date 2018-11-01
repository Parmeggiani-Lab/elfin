#!/usr/bin/env python3

#
# This script creates the atom model of an Elfin design solution
#   Input: a JSON file that describes the connectivity of a sol-
#       lution. 
#   Output: a CIF file define the positions of each atom of the 
#       input solution
#
# Version: v2
#
#

import glob
import argparse
import Bio.PDB
import Bio.SubsMat.MatrixInfo
import Bio.PDB.StructureBuilder

try:
    from utilities import *
    from pdb_utilities import *
except ImportError as e:
    from .utilities import *
    from .pdb_utilities import *


class Synthesiser:
    def __init__(
        self, 
        spec, 
        pdb_dir, 
        cappings_dir,
        metadata_dir, 
        show_fusion=False,
        disable_capping=False
    ):
        self.spec             = spec
        self.pdb_dir          = pdb_dir
        self.cr_dir           = cappings_dir
        self.show_fusion      = show_fusion
        self.disable_capping  = disable_capping

        self.si               = Bio.PDB.Superimposer()

        # Parse and convert capping repeat indicies into a dictionary
        self.capping_repeat_r_id_dict = {}
        for row in read_csv(metadata_dir + '/repeat_indicies.csv', delim=' '):
            self.capping_repeat_r_id_dict[row[0].split('.')[0].replace('DHR', 'D')] = [int(idx) for idx in row[1:]]

    def reset_residue_id(self):
        self.residue_id = 1

    def get_next_residue_id(self):
        rid = self.residue_id
        self.residue_id += 1
        return rid

    def get_capping(self, prim_res, cap_res, cr_r_ids, n_term=True):
        if self.disable_capping:
            return []

        n_cap_res = len(cap_res)

        # We could use as many shared residues as possible for alignment but
        # that could create gaps (jumps) due to some of the residues in
        # primary being affected by interface. Therefore here we use one eigth
        # like GenXDB does (repeat index range / 4 is 1/8 of the module).
        if n_term:
            match_start_idx = [i for (i,el) in enumerate(cap_res) if el.id[1] == cr_r_ids[0]][0]
            for cri in range(match_start_idx, n_cap_res):
                prim_idx = cri - match_start_idx
                if prim_res[prim_idx].resname != cap_res[cri].resname: 
                        break
                max_match_idx = cri

            max_align_len = int_floor((max_match_idx - match_start_idx) / 4)
            prim_align_res = prim_res[:max_align_len]
            cap_align_res = cap_res[match_start_idx:match_start_idx+max_align_len]
            real_cap_res = cap_res[:match_start_idx]

            print('N-Terminus capping align len: {}'.format(max_align_len))
        else:
            match_end_idx = [i for (i,el) in enumerate(cap_res) if el.id[1] == cr_r_ids[3]][0]
            for cri in reversed(range(match_end_idx + 1)):
                prim_idx = cri - match_end_idx - 1
                if prim_res[prim_idx].resname != cap_res[cri].resname: 
                        break
                min_match_idx = cri

            max_align_len = int_floor((match_end_idx - min_match_idx) / 4)
            prim_align_res = prim_res[-max_align_len:]
            cap_align_res = cap_res[match_end_idx-max_align_len+1:match_end_idx+1]
            real_cap_res = cap_res[match_end_idx+1:]

            print('C-Terminus capping align len: {}'.format(max_align_len))

        for i in range(len(prim_align_res)):
            assert(prim_align_res[i].resname == cap_align_res[i].resname)

        prim_atoms       = [al[0] for al in [[a for a in r.child_list if a.name == 'CA'] for r in prim_align_res]]
        cap_atoms        = [al[0] for al in [[a for a in r.child_list if a.name == 'CA'] for r in cap_align_res]]
        
        self.si.set_atoms(prim_atoms, cap_atoms)
        rot, tran = self.si.rotran

        for r in real_cap_res:
            r.transform(rot, tran)
        
        return real_cap_res

    def project_module_instance(
        self, 
        graph,
        chains, 
        node_a
    ):
        print('Processing node: id={}, name={}'.format(node_a['id'], node_a['name']))

        if not node_a['trim'][1]:
            # This is an end-node and end-node atoms are already covered by
            # their preceeding non-end node
            return chains

        single_a = read_pdb(self.pdb_dir + '/singles/' + node_a['name'] + '.pdb')
        resi_count_a = get_pdb_residue_count(single_a)

        cterm_nodes = cterm_nodes=[n for n in graph['nodes'] if n['id'] == node_a['cterm_node_id']]
        n_cterm_nodes = len(cterm_nodes)

        if n_cterm_nodes == 0:
            # end-nodes don't reach here - this only fires in cases where
            # there's only one single node in the entire chain
            print('Warning: \n'
                    '   Single-node chains still use transformation '
                    '   matricies calculated with double fusion. '
                    '   This might cause some inaccuracy when applied '
                    '   not on a trimmed double but a single module.')
            chain_id = graph['name'] + '_' + str(node_a['id'])
            residues = strip_residues(single_a)

        elif n_cterm_nodes == 1:
            node_b = cterm_nodes[0]

            singleB = read_pdb(self.pdb_dir + '/singles/' + node_b['name'] + '.pdb')
            double = read_pdb(self.pdb_dir + '/doubles/' + (node_a['name'] + '-' + node_b['name']) + '.pdb')

            resi_count_b = get_pdb_residue_count(singleB)
            resi_count_double = get_pdb_residue_count(double)

            chain_id = graph['name'] + '_' + str(node_a['id']) + '-' + str(node_b['id'])
            residues = strip_residues(double)

            # Compute trim residue indicies that minimise effect on 
            # atom position caused by double interfaces
            prefix_residues = []
            suffix_residues = []
            if node_a['trim'][0]:
                start_resi = int_floor(float(resi_count_a)/2)
            else:
                start_resi = 0

                # We expect an untrimmed end to be capped
                if node_a['cap'][0]:
                    cap_name = node_a['name'].split('_')[0]
                    cap_and_repeat = read_pdb(self.cr_dir + '/' + cap_name + '_NI.pdb')
                    prefix_residues = self.get_capping(
                        prim_res=residues, 
                        cap_res=strip_residues(cap_and_repeat), 
                        cr_r_ids=self.capping_repeat_r_id_dict[cap_name], 
                        n_term=True
                    )
                else:
                    print('Warning: untrimmed N-Terminus is not capped')
            
            if node_b['trim'][1]:
                end_resi = resi_count_double - int_ceil(float(resi_count_b)/2)
            else:
                end_resi = resi_count_double

                # We expect an untrimmed end to be capped
                if node_b['cap'][1]:
                    cap_name = node_b['name'].split('_')[-1]
                    cap_and_repeat = read_pdb(self.cr_dir + '/' + cap_name + '_IC.pdb')
                    suffix_residues = self.get_capping(
                        prim_res=residues, 
                        cap_res=strip_residues(cap_and_repeat), 
                        cr_r_ids=self.capping_repeat_r_id_dict[cap_name], 
                        n_term=False
                    )
                else:
                    print('Warning: untrimmed C-Terminus is not capped')

            residues = prefix_residues + residues[start_resi:end_resi] + suffix_residues
        else:
            raise ValueError('n_cterm_nodes = ' + str(n_cterm_nodes))

        # Steal residues from the PDB and put in our chain
        if self.show_fusion:
            curr_chain = Bio.PDB.Chain.Chain(chain_id)
            chains.append(curr_chain)
        else:
            curr_chain = chains[0]

        for r in residues:
            r.id = (r.id[0], self.get_next_residue_id(), r.id[2]) 
            r.transform(np.asarray(node_a['rot']), node_a['tran'])
            curr_chain.add(r)

    def make_chains(self, graph):
        chains = []

        if not self.show_fusion:
            chains.append(Bio.PDB.Chain.Chain(graph['name']))

        for node in graph['nodes']:
            self.project_module_instance(
                graph=graph,
                chains=chains, 
                node_a=node
            )

        return chains

    def run(self):
        model = Bio.PDB.Model.Model(0)
        si = Bio.PDB.Superimposer()

        if self.show_fusion:
            print('Note: show_fusion is on')

        self.reset_residue_id()
        for graph_idx in range(len(self.spec)):
            try:
                print('Processing graph: idx={}, name={}'.format(graph_idx, self.spec[graph_idx]['name']))
            except KeyError:
                print('Error!')
                print('Input file seems to be v1 format.')
                print('Use v1_design_convert.py to convert a v1 format file to v2 format.')
                exit(1)

            chains = self.make_chains(self.spec[graph_idx])
            for c in chains:
                model.add(c)

        sb = Bio.PDB.StructureBuilder.StructureBuilder()
        sb.init_structure('0')
        structure = sb.get_structure()
        structure.add(model)
        return structure

def parse_args(args):
    parser = argparse.ArgumentParser(description='Generate atom model in CIF format using output from Elfin core');
    parser.add_argument('spec_file')
    parser.add_argument('--out_file', default='')
    parser.add_argument('--pdb_dir', default='./resources/pdb_aligned/')
    parser.add_argument('--cappings_dir', default='./resources/pdb_relaxed/cappings')
    parser.add_argument('--metadata_dir', default='./resources/metadata/')
    parser.add_argument('--show_fusion', action='store_true')
    parser.add_argument('--disable_capping', action='store_true')
    return parser.parse_args(args)

def main(test_args=None):
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    spec_ext = args.spec_file[args.spec_file.rfind('.'):]

    if spec_ext == '.json':
        spec = read_json(args.spec_file)
    else:
        print('Unknown spec file type: {}'.format(spec_ext))
        exit()

    if len(spec) > 1:
        print('Warning: multi-chain feature is not well tested yet')

    struct = Synthesiser(
        spec, 
        args.pdb_dir,
        args.cappings_dir,
        args.metadata_dir,
        args.show_fusion,
        args.disable_capping
    ).run()

    if args.out_file == '':
        args.out_file = args.spec_file
    args.out_file = '.'.join(args.out_file.split('.')[:-1] + ['cif'])

    print('Saving to:', args.out_file)
    save_cif(struct=struct, path=args.out_file)

    # Todo: output coms for multiple chains
    #   Need to change csv format such that points
    # aren't connected between chains
    
    # coms = [node['tran'] for node in spec[0]['nodes']]
    # saveCsv(coms, args.out_file + '.csv')

if __name__ == '__main__':
    safe_exec(main)