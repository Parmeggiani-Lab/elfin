#!/usr/bin/env python3

#
# This script creates the CIF atom model from a design solution exported from
# elfin-ui in JSON format.
#

from collections import deque
from collections import namedtuple
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

def get_node(network, json_name):
    node = network[json_name]
    check_mod_type(node['module_type'])
    return node

TermIdentifierBase = namedtuple('TermIdentifierBase', ['ui_name', 'chain_id', 'term'])
class TermIdentifier(TermIdentifierBase):
    """Small class to hold terminus identifier data"""
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls, *args, **kwargs)
        check_term_type(self.term)
        return self
    
    def __repr__(self):
        return ':'.join((str(getattr(self, f)) for f in self._fields))

    def same_chain_as(self, other):
        return self.ui_name == other.ui_name and \
            self.chain_id == other.chain_id

ChainIdentifierBase = namedtuple('ChainIdentifierBase', ['src', 'dst'])
class ChainIdentifier(ChainIdentifierBase):
    """Small class to hold source and destination TermIdentifiers for a chain"""
    def __new__(cls, *args, **kwargs):
        self = super().__new__(cls, *args, **kwargs)
        assert self.src.term == 'n'
        assert self.dst.term == 'c'
        return self
    
    def __repr__(self):
        return '{}->{}'.format(self.src, self.dst)

# Returns a list of all leaf TermIdentifiers.
#
# A leaf is a node (whether single or hub) with only one occupied terminus.
def find_leaves(network, xdb):
    try:
        res = []

        for ui_name in network:
            node = get_node(network, ui_name)

            cl = node['c_linkage']
            nl = node['n_linkage']

            # A single can be an entry if it only has one busy terminus. When
            # entry is found, return the free (far end) terminus identifier.
            if len(cl) == 1 and not nl:
                res.append(TermIdentifier(ui_name, cl[0]['source_chain_id'], 'n'))
            elif len(nl) == 1 and not cl:
                res.append(TermIdentifier(ui_name, nl[0]['source_chain_id'], 'c'))
            elif not nl and not cl:
                mod_type = node['module_type']
                mod_name = node['module_name']
                chains = xdb['modules'][mod_type + 's'][mod_name]['chains']
                for c in chains:
                    if chains[c]['n']:
                        res.append(TermIdentifier(ui_name, c, 'n'))
                    if chains[c]['c']:
                        res.append(TermIdentifier(ui_name, c, 'c'))
            
        return res
    except KeyError as ke:
        print('KeyError:', ke)
        print('Probably bad input format.')
        exit()

# Walks the a chain starting with the src TermIdentifier, according to the
# network JSON object, yielding each TermIdentifier and next_linkage on the
# fly.
def walk_chain(network, src):
    ui_name, chain_id, term = src

    while True:
        node = get_node(network, ui_name)

        # Advance until either a hub or a single with dangling terminus is
        # encountered.
        next_linkages = [l for l in node[opposite_term(term) + '_linkage'] \
            if l['source_chain_id'] == chain_id]

        assert len(next_linkages) <= 1, \
            'Expected only next_linkages size <= 1 (since' \
            'each node has max 2 linkages, one N and one C).'

        term_iden = TermIdentifier(ui_name, chain_id, term)
        next_linkage = next_linkages[0] if next_linkages else None
        yield term_iden, next_linkage

        if not next_linkage:
            break

        ui_name, chain_id = \
            next_linkage['target_mod'], \
            next_linkage['target_chain_id']

# Walks the network and returns a generator of ChainIdentifiers.
# 
# This method guarantees that src->dst is in the direction of N->C.
def decompose_network(network, xdb, skip_unused=False):
    src_q = deque()
    visited = set()

    # Find entry node to begin walking the network with.
    leaves = find_leaves(network, xdb)
    assert leaves, 'No leave nodes for network.'

    src_q.extend(leaves)

    # Work until queue is empty.
    while src_q:
        src = src_q.popleft()
        if src in visited:
            # This could happen when termini identifiers on hubs are added
            # before the termini on the other end of those chains are popped
            # out of the queue.
            continue

        visited.add(src)
        mod_type = None

        chain_walker = walk_chain(network, src)
        for term_iden, next_linkage in chain_walker:
            ui_name, chain_id, term = term_iden
            node = get_node(network, ui_name)

            if not next_linkage:
                dst = TermIdentifier(ui_name, chain_id, opposite_term(term))
                if dst not in visited:
                    visited.add(dst)

                    srcdst = (src, dst) if term == 'n' else (dst, src)
                    chain_iden = ChainIdentifier(*srcdst)

                    yield chain_iden

                    if mod_type == 'hub':            
                        # Add unvisited components as new chain sources.
                        hub = xdb['modules']['hubs'][node['module_name']]
                        for hub_chain_id in hub['chains']:
                            hub_chain = hub['chains'][hub_chain_id]
                            for term in TERM_TYPES:
                                if hub_chain[term]: # If not dormant.
                                    iden = (ui_name, hub_chain_id, term)
                                    if iden not in visited:
                                        src_q.append(iden)
                break

            mod_type = node['module_type']
            if mod_type == 'hub':
                # This is a "bypass" hub, i.e. the current hub component has
                # interfaceable N and C terms, and the current chain goes
                # through it without ending here.
                #
                # In this case, check for unused components that might not
                # need to be placed since they aren't leaves nor connect to
                # any leaf nodes.
                hub = xdb['modules']['hubs'][node['module_name']]
                for hub_chain_id in hub['chains']:
                    if hub_chain_id == chain_id: continue
                    c_links = len([l for l in node['c_linkage'] \
                        if l['source_chain_id'] == hub_chain_id])
                    n_links = len([l for l in node['n_linkage'] \
                        if l['source_chain_id'] == hub_chain_id])

                    if c_links == n_links == 0:
                        if skip_unused:
                            print('Skipping unused hub component:', ui_name, hub_chain_id)
                        else:
                            srcdst = (TermIdentifier(ui_name, hub_chain_id, 'n'),
                                TermIdentifier(ui_name, hub_chain_id, 'c'))
                            yield ChainIdentifier(*srcdst)

def trim_terminus(residues, term):
    check_term_type(term)

    if term == 'n':
        print('TODO: Trim term', term)
    elif term == 'c':
        print('TODO: Trim term', term)

def cap_terminus(residues, term, network, xdb, term_iden):
    check_term_type(term)

    node = get_node(network, term_iden.ui_name)
    mod_type = node['module_type']

    if mod_type == 'single':
        if term == 'n':
            print('TODO: Cap single term', term)
        elif term == 'c':
            print('TODO: Cap single term', term)
    elif mod_type == 'hub':
        # I'm not sure how to cap hubs yet, because they have
        # different residue sequence vs. singles (and hence caps).
        #
        # If we were to cap hubs, we need to first check whether N
        # term is an open terminus in this hub.

        hub = xdb['modules']['hubs'][node['module_name']]
        chain = hub['chains'][term_iden.chain_id]
        if chain[term]:
            print('WARNING: capping for open hub interface is not supported ({})'.format(term_iden))
        else:
            # No need to cap a hub component term that is a closed interface.
            pass

ModInfo = namedtuple('ModInfo', ['mod_type', 'mod_name', 'res', 'res_n'])

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

    def deposit_chain(self, network, chain_iden):
        # n -src-> c ... n -dst-> c
        print('Deposit chain:', chain_iden)

        src, dst = chain_iden
        start_rid = self.residue_id
        chain = self.new_chain()

        chain_walker = walk_chain(network, src)
        ids = set()
        for term_iden, next_linkage in chain_walker:
            print('Deposit {}->{}'.format(repr(term_iden), next_linkage['target_mod'] \
                if next_linkage else None))

            a_ui_name, a_chain_id, _ = term_iden  # term is always 'n'
            a_node = get_node(network, a_ui_name)
            a_info = self.get_mod_info(a_node, a_chain_id)
            a_info2 = self.get_mod_info(a_node, a_chain_id)

            residues = a_info.res
            # if not chain and not next_linkage:
            #     # The chain consists of just one module.
            #     residues = a_info.res
            # elif next_linkage:
            #     b_ui_name, b_chain_id = next_linkage['target_mod'], \
            #         next_linkage['target_chain_id']
            #     b_node = get_node(network, b_ui_name)
            #     b_info = self.get_mod_info(b_node, b_chain_id)

            #     types = (a_info.mod_type, b_info.mod_type) 

            #     if types == ('single', 'single'):
            #         dbl_name = a_info.mod_name + '-' + b_info.mod_name
            #         dbl_pdb = read_pdb(self.pdb_dir + '/doubles/' + dbl_name + '.pdb')
            #         dbl_res = copy_residues(dbl_pdb)
            #         dbl_res_n = len(dbl_res)

            #         residues = [r.copy() for r in dbl_res]

            #     elif types == ('hub', 'single'):
            #         print('TODO: hub-single deposit')
            #     elif types == ('single', 'hub'):
            #         print('TODO: single-hub deposit')
            #     else:
            #         raise ValueError('Unknown type tuple:', types)

            # # Compute trim residue indicies that minimise effect on atom
            # # position caused by double interfaces
            # prefix_residues = []
            # suffix_residues = []
            # if chain:
            #     # Midway through the chain - always trim N term.
            #     start_resi = float(a_info.res_n)//2
            # else:
            #     # Start of chain on the N side - cap N term.  
            #     start_resi = 0

            #     cap_name = a_info.mod_name.split('_')[0]
            #     cap_and_repeat = read_pdb(self.cr_dir + '/' + cap_name + '_NI.pdb')
            #     prefix_residues = self.get_capping(
            #         prim_res=residues, 
            #         cap_res=get_residues(cap_and_repeat), 
            #         cr_r_ids=self.capping_repeat_idx[cap_name], 
            #         term='n'
            #     )
            
            # if node_b['trim']['c']:
            #     end_resi = dbl_res_n - int_ceil(float(b_info.res_n)/2)
            # else:
            #     end_resi = dbl_res_n

            #     # We expect an untrimmed end to be capped
            #     if node_b['cap'][1]:
            #         cap_name = node_b['name'].split('_')[-1]
            #         cap_and_repeat = read_pdb(self.cr_dir + '/' + cap_name + '_IC.pdb')
            #         suffix_residues = self.get_capping(
            #             prim_res=residues, 
            #             cap_res=get_residues(cap_and_repeat), 
            #             cr_r_ids=self.capping_repeat_idx[cap_name], 
            #             term='c'
            #         )
            #     else:
            #         print('Warning: untrimmed C-Terminus is not capped')

            # residues = prefix_residues + residues[start_resi:end_resi] + suffix_residues

            # # Cap or trim N term?
            # if chain:
            #     # Midway through the chain - always trim N term.
            #     trim_terminus(residues, 'n')
            # else:
            #     # Start of chain on the N side - cap N term. 
            #     cap_terminus(residues, 'n', network, self.xdb, term_iden)

            # # Cap or trim C term?
            # if next_linkage:
            #     # There's a next node - always trim C term.
            #     trim_terminus(residues, 'c')
            # else:
            #     # There's no next node - cap C term.
            #     cap_terminus(residues, 'c', network, self.xdb, term_iden)

            for r in residues:
                r.id = (r.id[0], self.next_residue_id(), r.id[2])
                rot = np.transpose(np.asarray(a_node['rot']))
                r.transform(rot, a_node['tran'])
                chain.add(r)

            if self.show_fusion:
                # curr_chain = Bio.PDB.Chain.Chain(chain_id)
                # chains.append(curr_chain)
                print('TODO: show_fusion')

        print('')
        self.model.add(chain)

    def get_capping(self, prim_res, cap_res, cr_r_ids, term):
        if self.disable_capping:
            return []

        check_term_type(term)

        cap_res_n = len(cap_res)
        prim_res_n = len(prim_res)
        prim_align_res = []

        if term == 'n':
            # <Capping Residues> <Match Start ... End> <Rest of Primary...>

            # Find residue index at which the residue id[1] matches capping
            # start index. Residue id often does not start from 1 and is never
            # 0-based.
            min_match_idx = [i for (i,el) in enumerate(cap_res) \
                if el.id[1] == cr_r_ids[0]][0]

            # Find max match index, which is either the number of capping
            # residues or the first index at which residue sequences diverge.
            id_range = range(min_match_idx, cap_res_n)
            for i in id_range:
                prim_r = prim_res[i]
                if prim_r.resname != cap_res[i].resname: 
                        break
                max_match_idx = i
                prim_align_res.append(prim_r)

            real_cap_res = cap_res[:min_match_idx]

        elif term == 'c':
            # <Rest of Primary...> <Match Start ... End> <Capping Residues>
            max_match_idx = [i for (i,el) in enumerate(cap_res) \
                if el.id[1] == cr_r_ids[3]][0]

            id_range = range(0, max_match_idx + 1)
            for i in reversed(id_range):
                prim_id = prim_res_n - 1 - (max_match_idx - i)
                prim_r = prim_res[prim_id]
                if prim_r.resname != cap_res[i].resname: 
                        break
                min_match_idx = i
                prim_align_res.append(prim_r)

            prim_align_res.reverse()
            real_cap_res = cap_res[max_match_idx+1:]

        cap_align_res = cap_res[min_match_idx:max_match_idx+1]
        align_len = (max_match_idx - min_match_idx) // 4
        print('------DEBUG: {} term capping align len: {}'.format(term, align_len))

        prim_atoms = [r['CA'] for r in prim_align_res]
        cap_atoms = [r['CA'] for r in cap_align_res]
        
        self.si.set_atoms(prim_atoms, cap_atoms)
        rot, tran = self.si.rotran

        result = []
        for r in real_cap_res:
            rr = r.copy()
            rr.transform(rot, tran)
            result.append(rr)
        
        return result

    def get_mod_info(self, node, chain_id):
        mod_type = node['module_type']
        mod_name = node['module_name']

        # Obtain module residues
        pdb = read_pdb(self.pdb_dir + '/' + mod_type + 's/' + mod_name + '.pdb')
        res = copy_residues(pdb, chain_id)
        res_n = len(res)

        return ModInfo(mod_type, mod_name, res, res_n)

    def deposit_chains(self, network):
        chain_iden_gen = decompose_network(network, self.xdb, self.skip_unused)
        for chain_iden in chain_iden_gen:
            self.deposit_chain(network, chain_iden)

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

        pause_code()
        return structure

if __name__ == '__main__':
    safe_exec(main)