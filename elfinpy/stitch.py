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
    desc = 'Create CIF atom model from design solution JSON exported by elfin-ui.'
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

        struct = Stitcher(
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
# A leaf is a terminus that is either unoccupied or on a hub node.
def find_leaves(network, xdb):
    try:
        res = []

        for ui_name in network:
            node = get_node(network, ui_name)

            mod_type = node['module_type']
            mod_name = node['module_name']
            chains = xdb['modules'][mod_type + 's'][mod_name]['chains']
            cl = node['c_linkage']
            nl = node['n_linkage']

            if mod_type == 'hub':
                for c in chains:
                    if chains[c]['n']:
                        res.append(TermIdentifier(ui_name, c, 'c'))
                    if chains[c]['c']:
                        res.append(TermIdentifier(ui_name, c, 'n'))
            else:  # Guaranteed to be 'single' thanks to get_node()
                if not nl:
                    res.append(TermIdentifier(ui_name, cl[0]['source_chain_id'], 'n'))
                if not cl:
                    res.append(TermIdentifier(ui_name, nl[0]['source_chain_id'], 'c'))
            
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

ModInfo = namedtuple('ModInfo', ['mod_type', 'mod_name', 'res', 'res_n'])

def transform_residues(res, rot, tran):
    for r in res:
        for a in r:
            # Do the transform manually because BioPython has non
            # standard multiplication order.
            a.coord = rot.dot(a.coord) + tran

# Blend residue lists M = (1-w)M + wF, where M is an atom coordinate in
# moving_res, F is an atom coordinate in fixed_res, and w is the corresponding
# weight in weights.
#
# Also removes dirty atoms. If residues are not the same (name), only backbone
# atoms are blended.
def blend_residues(moving_res, fixed_res, weights):
    assert len(moving_res) == len(fixed_res)
    assert len(moving_res) == len(weights)

    for m, f, w in zip(moving_res, fixed_res, weights):
        # Remove dirty atoms. They seem crop up in the process of
        # optimizing PDBs even if preprocess.py already removed them once.
        #
        # Also remove atoms not in fixed residue - this is only known to
        # happen to CYS (HG) and HIS (HE1/HE2).
        to_remove = [a for a in m if a.name in DIRTY_ATOMS or a.name not in f]
        for da in to_remove:
            m.detach_child(da.name)

        # sidechain_atoms = [a for a in m if a.name not in f]
        # for sa in sidechain_atoms:
        #     m.detach_child(sa.name)

        # Compute new position based on combination of two positions.
        compute_coord = lambda a, b : (1-w)*a.coord + w*b.coord

        for ma in m:
            if m.resname == f.resname:
                assert ma.name in f  # Identical residues should have the same atoms
                ma.coord = compute_coord(ma, f[ma.name])
            else:
                # Only modify backbone atoms.
                if ma.name in BACKBONE_NAMES and \
                    ma.name in f:
                    ma.coord = compute_coord(ma, f[ma.name])

class Stitcher:
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
        atom_chain = self.new_chain()

        # Build context to pass to subroutines.
        context = lambda: 0
        context.atom_chain = atom_chain
        context.network = network
        context.last_node = None

        chain_walker = walk_chain(network, src)
        for term_iden, next_linkage in chain_walker:
            context.term_iden = term_iden
            context.next_linkage = next_linkage

            print('Deposit {}->{}'.format(repr(term_iden), next_linkage['target_mod'] \
                if next_linkage else None))

            context.node = get_node(network, term_iden.ui_name)
            context.mod_info = self.get_mod_info(context.node, term_iden.chain_id)

            context.pref_res = []
            context.main_res = [r.copy() for r in context.mod_info.res]
            context.suff_res = []

            if context.last_node:
                # Midway through the chain - always displace N term.
                self.displace_terminus(context, 'n')
            else:
                # Start of chain on the N side - cap N term. 
                self.cap_terminus(context, 'n')

            if next_linkage:
                # There's a next node - always displace C term.
                self.displace_terminus(context, 'c')
            else:
                # There's no next node - cap C term.
                self.cap_terminus(context, 'c')

            all_res = context.pref_res + context.main_res + context.suff_res
            for r in all_res:
                r.id = (r.id[0], self.next_residue_id(), r.id[2])
                rot = np.transpose(np.asarray(context.node['rot']))
                r.transform(rot, context.node['tran'])
                atom_chain.add(r)

            if self.show_fusion:
                # curr_chain = Bio.PDB.Chain.Chain(chain_id)
                # chains.append(curr_chain)
                print('TODO: show_fusion')

            context.last_node = context.node
            context.last_term_iden = term_iden

        print('')
        self.model.add(atom_chain)

    def cap_terminus(self, deposit_context, term):
        check_term_type(term)

        # Unpack context.
        mod_info = deposit_context.mod_info
        residues = deposit_context.main_res
        chain_id = deposit_context.term_iden.chain_id

        if mod_info.mod_type == 'single':
            cap_name = mod_info.mod_name.split('_')[-1]
            print('Cap({}): {}'.format(term, cap_name))
            pdb_path = '{}/{}_{}.pdb'.format(self.cr_dir, cap_name,
            'NI' if term == 'n' else 'IC')
            cap_and_repeat = read_pdb(pdb_path)

            cap_res = self.get_capping(
                prime_res=residues, 
                cap_res=get_residues(cap_and_repeat), 
                cr_r_ids=self.capping_repeat_idx[cap_name], 
                term=term
            )

            if term == 'n':
                deposit_context.pref_res = cap_res
            else:
                deposit_context.suff_res
        elif mod_info.mod_type == 'hub':
            # If we were to cap hubs, we need to first check whether N
            # term is an open terminus in this hub.
            hub = self.xdb['modules']['hubs'][mod_info.mod_name]
            chain = hub['chains'][chain_id]
            if chain[term]:
                print('TODO: Cap hub term', term)
            else:
                # No need to cap a hub component term that is a closed interface.
                pass

    # Computes the capping residues. Displaces primary residues (thus modifies
    # the prime_res parameter).
    def get_capping(self, prime_res, cap_res, cr_r_ids, term):
        if self.disable_capping:
            return []

        check_term_type(term)

        cap_res_n = len(cap_res)
        prim_res_n = len(prime_res)

        # Find residue index at which the residue id[1] matches capping
        # start index. Residue id often does not start from 1 and is never
        # 0-based.
        rid_range = tuple(cr_r_ids[:2]) if term == 'n' else tuple(cr_r_ids[2:])

        for i, el in enumerate(cap_res):
            if el.id[1] == rid_range[0]:
                match_start = i
                break
        else:
            raise ValueError('Could not find residue index {}'.format(
                rid_range[0]))

        match_len = rid_range[1] - rid_range[0]
        match_end = match_start + match_len

        # N: match left, C: match right
        prime_align_res = prime_res[:match_len] if term =='n' else prime_res[-match_len:]
        cap_align_res = cap_res[match_start:match_end]

        prim_atoms = [r['CA'] for r in prime_align_res]
        cap_atoms = [r['CA'] for r in cap_align_res]
        
        self.si.set_atoms(prim_atoms, cap_atoms)
        rot, tran = self.si.rotran

        result = []
        cap_protrude_res = cap_res[:match_start] + cap_res[match_end+1:]
        for r in cap_protrude_res:
            rr = r.copy()
            rr.transform(rot, tran)
            result.append(rr)

        # Also transform cap align res to the right frame.
        for r in cap_align_res:
            r.transform(rot, tran)

        # Displace prime_res using linear weights, the same method as
        # displace_termins().

        # Linear weights (0, 1].
        disp_w = [i/match_len for i in range(1, match_len + 1)]
        if term == 'n':
            disp_w.reverse()  # Want [1, 0) for N term.

        blend_residues(prime_align_res, cap_align_res, disp_w)

        return result

    def displace_terminus(self, deposit_context, term):
        check_term_type(term)

        if term == 'n':
            assert deposit_context.last_node

            # Node A is on the C end, so we get the N end node, and swap
            # order.
            a_node = deposit_context.last_node
            a_chain_id = deposit_context.last_term_iden.chain_id
            a_info = self.get_mod_info(a_node, a_chain_id)

            b_node = deposit_context.node
            b_chain_id = deposit_context.term_iden.chain_id
            b_info = deposit_context.mod_info

        elif term == 'c':
            next_linkage = deposit_context.next_linkage
            assert next_linkage

            # Node A is on the N end, so we get the C end node.
            a_node = deposit_context.node
            a_chain_id = deposit_context.term_iden.chain_id
            a_info = deposit_context.mod_info

            b_ui_name, b_chain_id = next_linkage['target_mod'], \
                next_linkage['target_chain_id']
            b_node = get_node(deposit_context.network, b_ui_name)
            b_info = self.get_mod_info(b_node, b_chain_id)

        types = (a_info.mod_type, b_info.mod_type)
        if types == ('single', 'single'):
            a_single_name = a_info.mod_name
            b_single_name = b_info.mod_name

        elif types == ('hub', 'single'):
            hub = self.xdb['modules']['hubs'][a_info.mod_name]
            a_single_name = hub['chains'][a_chain_id]['single_name']
            b_single_name = b_info.mod_name

        elif types == ('single', 'hub'):
            a_single_name = a_info.mod_name
            hub = self.xdb['modules']['hubs'][b_info.mod_name]
            b_single_name = hub['chains'][b_chain_id]['single_name']
        else:
            raise ValueError('Unknown type tuple:', types)

        a_single_len = self.get_single_len(a_single_name)

        dbl_name = a_single_name + '-' + b_single_name
        dbl_pdb = read_pdb(self.pdb_dir + '/doubles/' + dbl_name + '.pdb')
        dbl_res = get_residues(dbl_pdb)

        main_res = deposit_context.main_res
        disp_n = len(dbl_res) // 2

        # Linear weights (0, 1].
        disp_w = [i/disp_n for i in range(1, disp_n + 1)]
        
        if term == 'n':
            if b_info.mod_type == 'hub':
                # Lift double (in A frame) to hub arm frame with A at the
                # arm's tip.
                tx_id = self.xdb['modules']['singles'][a_single_name] \
                    ['chains'][a_chain_id]['c'][b_info.mod_name][b_chain_id]
                tx = self.xdb['n_to_c_tx'][tx_id]
                rot = np.asarray(tx['rot'])
                tran = np.asarray(tx['tran'])
            else:  # Guaranteed to be 'single' thanks to get_node()
                # Drop double to B frame.
                rot, tran = self.get_drop_tx(a_single_name, b_single_name)
            
            transform_residues(dbl_res, rot, tran)
            
            # Displace N term residues (first half of main_res) based on
            # linear weights. In the double, start from B module.
            #
            # main_res:                  [n ... | ... c]
            # disp_w:                    [1....0]
            # dbl:      [n ... | ... c]  [n ... | ... c]
            main_disp = main_res[:disp_n]
            dbl_part = dbl_res[a_single_len:a_single_len+disp_n]
            disp_w.reverse()  # Make it 1 -> 0
        elif term == 'c':            
            # Displace C term residues (second half of main_res) based on
            # linear weights. In the double, start from end of A module and go
            # backwards.
            #
            # main_res: [n ... | ... c]
            # disp_w:          [0....1]
            # dbl:      [n ... | ... c]  [n ... | ... c]
            if a_info.mod_type == 'hub':
                # Step 1: Drop double to B frame.
                rot, tran = self.get_drop_tx(a_single_name, b_single_name)

                # Step 2: Lift double (in B frame) to hub arm frame.
                tx_id = self.xdb['modules']['hubs'][a_info.mod_name] \
                    ['chains'][a_chain_id]['c'][b_single_name][b_chain_id]
                tx = self.xdb['n_to_c_tx'][tx_id]
                hub_rot = np.asarray(tx['rot'])
                hub_tran = np.asarray(tx['tran'])

                tran = hub_rot.dot(tran) + hub_tran
                rot = hub_rot.dot(rot)
                transform_residues(dbl_res, rot, tran)

            main_disp = main_res[-disp_n:]
            dbl_part = dbl_res[a_single_len-disp_n:a_single_len]

        blend_residues(main_disp, dbl_part, disp_w)

    def get_drop_tx(self, a_single_name, b_single_name):
        a_chains = self.xdb['modules']['singles'][a_single_name]['chains']

        assert len(a_chains) == 1
        a_chain_id = list(a_chains.keys())[0]
        a_b_chains = a_chains[a_chain_id]['c'][b_single_name]

        assert len(a_b_chains) == 1
        b_chain_id = list(a_b_chains.keys())[0]
        tx_id = a_b_chains[b_chain_id]
       
        tx = self.xdb['n_to_c_tx'][tx_id]

        # Inverse tx because dbgen.py computes the tx that takes the
        # single B module to part B inside double.
        rot = np.transpose(tx['rot'])
        tran = rot.dot(-np.asarray(tx['tran']))

        return rot, tran

    # Returns the number of residues in a single module.
    def get_single_len(self, mod_name):
        chains = self.xdb['modules']['singles'][mod_name]['chains']
        assert len(chains) == 1

        chain_id = list(chains.keys())[0]
        return chains[chain_id]['n_residues']

    def get_mod_info(self, node, chain_id):
        mod_type = node['module_type']
        mod_name = node['module_name']

        # Obtain module residues.
        pdb = read_pdb(self.pdb_dir + '/' + mod_type + 's/' + mod_name + '.pdb')
        res = get_residues(pdb, chain_id)
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

        return structure

if __name__ == '__main__':
    safe_exec(main)