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
    parser.add_argument('--xdb', default='./resources/xdb.json')
    parser.add_argument('--pdb_dir', default='./resources/pdb_aligned/')
    parser.add_argument('--cappings_dir', default='./resources/pdb_relaxed/cappings')
    parser.add_argument('--metadata_dir', default='./resources/metadata/')
    parser.add_argument('--show_fusion', action='store_true')
    parser.add_argument('--disable_capping', action='store_true')
    parser.add_argument('--skip_unused', action='store_true')
    return parser.parse_args(args)

def main(test_args=None):
    args = parse_args(sys.argv[1:] if test_args is None else test_args)

    input_ext = args.input_file[args.input_file.rfind('.'):].lower()

    if input_ext == '.json':
        spec = read_json(args.input_file)
        xdb = read_json(args.xdb)

        struct = Depositor(
            spec, 
            xdb,
            args.pdb_dir,
            args.cappings_dir,
            args.metadata_dir,
            args.show_fusion,
            args.disable_capping,
            args.skip_unused
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
        if not spec['networks']:
            return 'Spec file has no module networks.'

        if 'pg_networks' in spec:
            n_pgn = len(spec['pg_networks'])
            if n_pgn > 0:
                return \
                    'Spec file has {} path guide networks. It should have zero.'\
                    .format(n_pgn)

def find_leaves(network):
    # Find all leaf terminus identifier. A leaf node is one that has only one
    # occupied terminus.

    try:
        res = []

        for ui_name in network:
            uimod = network[ui_name]
            mod_type = uimod['module_type']

            cl = uimod['c_linkage']
            nl = uimod['n_linkage']

            # A single can be an entry if it only has one busy terminus.
            # When entry is found, return the outward terminus identifier.
            if len(cl) == 1 and not nl:
                res.append((ui_name, cl[0]['source_chain_id'], 'c'))
            elif len(nl) == 1 and not cl:
                res.append((ui_name, nl[0]['source_chain_id'], 'n'))
            
        return res
    except KeyError as ke:
        print('KeyError:', ke)
        print('Probably bad input format.')
        exit()

class Depositor:
    def __init__(
        self, 
        spec,
        xdb,
        pdb_dir, 
        cappings_dir,
        metadata_dir, 
        show_fusion=False,
        disable_capping=False,
        skip_unused=False,
    ):
        spec_complaint = validate_spec(spec)
        if spec_complaint:
            print('Error:', spec_complaint)
            exit()

        self.spec             = spec
        self.xdb              = xdb
        self.pdb_dir          = pdb_dir
        self.cr_dir           = cappings_dir
        self.show_fusion      = show_fusion
        self.disable_capping  = disable_capping
        self.skip_unused      = skip_unused
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
        if not chain:
            if module['module_type'] == 'single':
                # Cap on the unoccupied terminus.
                ...

            # Trim on occupied teriminus.
            ...

            # If module is a hub, then?
        else:
            # Trim on incoming terminus (chain_info).
            ...

            if module['module_type'] == 'single':
                if not module['c_linkage']:
                    # Cap on c terminus.
                    ...
                elif not module['n_linkage']:
                    # Cap on n terminus.
                    ...

            # If module is a hub, then?


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

    def deposit_chain(self, network, src, dst):
        # n -src-> c ... n -dst-> c
        assert(src and src[2] == 'c')
        assert(dst and dst[2] == 'n')

        print('Deposit:', src, dst)
        start_rid = self.residue_id
        chain = self.new_chain()


        self.model.add(chain)

    def deposit_chains(self, network):
        for chain_iden in self.decompose_network(network):
            self.deposit_chain(network, *chain_iden)

    def decompose_network(self, network):
        # Walks the network and returns a list of chain identifiers (src, dst),
        # where src and dst are outgoing and incoming terminus identifiers
        # respectively (ui_name, chain_id, term_type). This method guarantees that
        # src->dst is in the direction of N->C.
        #
        # A special case arises with unoccupied hub components or single modules.
        # These will have either src and dst being the same.

        src_q = deque()
        visited = set()

        # Find entry node to begin walking the network with.
        leaves = find_leaves(network)
        if not leaves:
            raise ValueError('No leave nodes for network.')

        src_q.extend(leaves)

        for src in src_q:
            print('Leaves:', src)

        # Work until queue is empty.
        res = []
        while src_q:
            src = src_q.popleft()
            if src in visited:
                # This could happen when termini identifiers on hubs are added
                # before the termini on the other end of those chains are
                # popped out of the queue.
                continue
            visited.add(src)

            ui_name, chain_id, term = src

            while True:
                node = network[ui_name]

                # Advance until either a hub or a single with dangling terminus is
                # encountered.
                next_linkage = [l for l in node[term + '_linkage'] \
                    if l['source_chain_id'] == chain_id]
                if not next_linkage:
                    break

                mod_type = node['module_type']
                if mod_type == 'hub':
                    # This is a bypass hub. Check for unused chains, which
                    # might not need to be placed since they aren't leaves nor
                    # do they connect to any leaf nodes.
                    hub = self.xdb['modules']['hubs'][node['module_name']]
                    for hub_chain_id in hub['chains']:
                        if hub_chain_id == chain_id: continue
                        c_links = len([l for l in node['c_linkage'] \
                            if l['source_chain_id'] == hub_chain_id])
                        n_links = len([l for l in node['n_linkage'] \
                            if l['source_chain_id'] == hub_chain_id])

                        if c_links == n_links == 0:
                            if self.skip_unused:
                                print('Skipping unused chain:', ui_name, hub_chain_id)
                            else:
                                res.append(((ui_name, hub_chain_id, 'c'),
                                    (ui_name, hub_chain_id, 'n')))

                assert(len(next_linkage) == 1)
                ui_name, chain_id = \
                    next_linkage[0]['target_mod'], \
                    next_linkage[0]['target_chain_id']

            dst = (ui_name, chain_id, opposite_term(term))
            if dst not in visited:
                visited.add(dst)
                res.append((src, dst) if term == 'c' else (dst, src))

                if mod_type == 'hub':            
                    # Add unvisited components as new chain sources.
                    hub = self.xdb['modules']['hubs'][node['module_name']]
                    for hub_chain_id in hub['chains']:
                        hub_chain = hub['chains'][hub_chain_id]
                        for term in TERM_TYPES:
                            if hub_chain[term]: # If not dormant.
                                iden = (ui_name, hub_chain_id, term)
                                if iden not in visited:
                                    src_q.append(iden)

        return res

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