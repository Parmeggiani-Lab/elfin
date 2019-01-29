#!/usr/bin/env python3

#
# This script creates the CIF atom model from a design solution exported from
# elfin-ui in JSON format.
#

from collections import deque
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


def parse_args(args):
    desc = 'Deposit CIF atom model from design solution JSON exported by elfin-ui.'
    parser = argparse.ArgumentParser(description=desc);
    parser.add_argument('input_file')
    parser.add_argument('--out_file', default='')
    parser.add_argument('--pdb_dir', default='./resources/pdb_aligned/')
    parser.add_argument('--cappings_dir', default='./resources/pdb_relaxed/cappings')
    parser.add_argument('--metadata_dir', default='./resources/metadata/')
    parser.add_argument('--show_fusion', action='store_true')
    parser.add_argument('--disable_capping', action='store_true')
    return parser.parse_args(args)

def main(test_args=None):
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    input_ext = args.input_file[args.input_file.rfind('.'):].lower()

    if input_ext == '.json':
        spec = read_json(args.input_file)
        struct = Depositor(
            spec, 
            args.pdb_dir,
            args.cappings_dir,
            args.metadata_dir,
            args.show_fusion,
            args.disable_capping
        ).run()

        if args.out_file == '':
            args.out_file = args.input_file
        args.out_file = '.'.join(args.out_file.split('.')[:-1] + ['cif'])

        print('Saving to:', args.out_file)
        save_cif(struct=struct, path=args.out_file)
    else:
        print('Unknown input file type: \"{}\"'.format(input_ext))
        exit()

def validate_spec(spec):
    if 'networks' not in spec:
        return 'No networks object in spec.'
    else:
        if len(spec['networks']) == 0:
            return 'Spec file has no module networks.'

        if 'pg_networks' in spec:
            n_pgn = len(spec['pg_networks'])
            if n_pgn > 0:
                return \
                    'Spec file has {} path guide networks. It should have zero.'\
                    .format(n_pgn)

def extract_linkage(linkage):
    return [(l['source_chain_id'], l['terminus'].lower()) for l in linkage]

class Depositor:
    def __init__(
        self, 
        spec, 
        pdb_dir, 
        cappings_dir,
        metadata_dir, 
        show_fusion=False,
        disable_capping=False
    ):
        spec_complaint = validate_spec(spec)
        if spec_complaint:
            print('Error:', spec_complaint)
            exit()

        self.spec             = spec
        self.pdb_dir          = pdb_dir
        self.cr_dir           = cappings_dir
        self.show_fusion      = show_fusion
        self.disable_capping  = disable_capping
        self.si               = Bio.PDB.Superimposer()
        self.chain_id         = 0

        # Parse and convert capping repeat indicies into a dictionary
        self.capping_repeat_idx = {}
        meta_csv = read_csv(metadata_dir + '/repeat_indicies.csv', delim=' ')
        for row in meta_csv:
            mod_name = row[0].split('.')[0].replace('DHR', 'D')
            self.capping_repeat_idx[mod_name] = \
                [int(idx) for idx in row[1:]]

    def place_module(self, chain, module, chain_info):
        print('Processing node: id={}, name={}'.format(node_a['id'], node_a['name']))

        if not node_a['trim']['c']:
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
            if node_a['trim']['n']:
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
            
            if node_b['trim']['c']:
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

    def deposit_chain(self, network, src_name, src_chain_info):
        start_rid = self.residue_id

        chain = new_chain()

        # Chain info is (Chain Name, Terminus, Trim)
        next_mod = src_name, src_chain_info

        while True:
            next_name, next_chain_info = next_mod
            next_mod = self.place_module(chain, network[next_name], next_chain_info)

            if not next_mod:
                raise ValueError('TODO: find chain info to return')
                # If next next_name is a hub, inactive components need to be
                # placed and capped.
                break

        _, _, src_cap = src_chain_info
        if src_term == 'c':
            # Need to reverse residue ids because N-to-C residue ids should be
            # in ascending order.
            raise ValueError('TODO: reverse reside ids')

        self.model.add(chain)

    def find_entry(self, network):
        # Find an arbitrary entry node and return the terminus direction
        # leading to subsequent nodes.

        try:
            for ui_name in network:
                uimod = network[ui_name]
                mod_type = uimod['module_type']
                    
                cl = uimod['c_linkage']
                nl = uimod['n_linkage']

                if mod_type == 'single':
                    # A single can be an entry if it only has one busy interface.
                    if len(cl) == 0 and len(nl) == 0:
                        raise ValueError('Orphant module: \"{}\"'.format(ui_name))

                    if len(cl) == 0:
                        assert(len(nl) == 1)
                        return ui_name, (nl[0]['source_chain_id'], nl[0]['terminus'].lower())
                    elif len(nl) == 0:
                        assert(len(cl) == 1)
                        return ui_name, (cl[0]['source_chain_id'], cl[0]['terminus'].lower())

        except KeyError as ke:
            print('KeyError:', ke)
            print('Probably bad input format.')
            exit()

    def deposit_chains(self, network):
        src_q = deque()
        visited = set()

        # Find entry node to begin walking the network with.
        entry = self.find_entry(network)
        if not entry:
            return 'Could not find entry node for network.'

        src_q.append(*entry)

        # Work until queue is empty.
        while not len(src_q) == 0:
            src = src_q.popleft()
            visited.add(src)

            dst_name, dst_chain_info = self.deposit_chain(network, *src)

            # Add unvisited branch as sources.
            for dt in dst_chain_info:
                if (dst_name, dst_chain_info) not in visited:
                    src_q.append((dst_name, dt))

    def new_chain(self):
        return Bio.PDB.Chain.Chain(self.next_chain_id())

    def next_chain_id(self):
        cid = str(self.chain_id)
        self.chain_id += 1
        return cid

    def reset_residue_id(self):
        self.residue_id = 1

    def next_residue_id(self):
        rid = self.residue_id
        self.residue_id += 1
        return rid

    def run(self):
        self.reset_residue_id()
        self.model = Bio.PDB.Model.Model(0)

        if self.show_fusion:
            print('Note: show_fusion is on')

        networks = self.spec['networks']
        for nw_name in networks:
            print('Processing network \"{}\"'.format(nw_name))
            complaint = self.deposit_chains(networks[nw_name])
            if complaint:
                print('Error: {}', complaint)
                exit()

        # Create output
        sb = Bio.PDB.StructureBuilder.StructureBuilder()
        sb.init_structure('0')
        structure = sb.get_structure()
        structure.add(self.model)
        return structure

if __name__ == '__main__':
    safe_exec(main)