import Bio.PDB

def get_pdb_residue_count(pdb):
    '''Returns the residue count of a Bio.PDB.Structure.Structure.'''
    return sum([len(c.child_list) for c in pdb.child_list[0].child_list])

def get_chain_residue_count(chain):
    '''Returns the residue count of a Bio.PDB.Structure.Structure.'''
    return len(chain.child_list)

def strip_residues(pdb, chain_ids=None):
    '''
    Returns returns residues removed from a PDB.

    Args:
    - pdb - Bio.PDB.Structure.Structure.
    - chain_ids - strip residues from these specific chain_ids only.

    Returns:
    - residues - a list of Bio.PDB.Residue.Residue.
    '''
    residues = []
    for model in pdb:
        for chain in model:
            if chain_ids == None or chain.id in chain_ids:
                tmp = list(chain.child_list)
                for r in tmp:
                    chain.detach_child(r.id)
                residues += tmp
    return residues

def get_chain(struct, chain_id='A'):
    '''Returns a specific chain from a Bio.PDB.Structure.Structure.'''
    return struct.child_list[0].child_dict[chain_id]

def get_chains(struct):
    '''Returns all chains of a Bio.PDB.Structure.Structure.'''
    return struct.child_list[0].child_list;

def read_pdb(
        read_path,
        pdb_name=None
):
    '''
    Reads a PDB file.

    Args:
    - read_path - PDB string file path to read from.
    - pdb_name - a string to set as the name of the Bio.PDB.Structure.Structure.
    
    Returns:
    - structure - Bio.PDB.Structure.Structure.
    '''
    if pdb_name == None:
        pdb_name = read_path.split('/')[-1].replace('.', '_')
    parser = Bio.PDB.PDBParser(PERMISSIVE=False)
    structure = parser.get_structure(pdb_name, read_path)
    return structure

def save_cif(**kwargs):
    '''
    Saves a Bio.PDB.Structure.Structure as a CIF file. Does not automatically
    append .cif extension.

    Args:
    - struct - Bio.PDB.Structure.Structure to be saved.
    - save_path - CIF string file path.
    '''
    struct = kwargs.pop('struct')
    save_path = kwargs.pop('save_path')

    with open(save_path, 'w') as file:
        io = Bio.PDB.mmcifio.MMCIFIO()
        io.set_structure(struct)
        io.save(file)
        # Temporary fix for CIF files not getting parsed properly by Rosetta: add
        # a dummy section at the end. ("Note that the final table in the cif file
        # may not be recognized - adding a dummy entry (like `_citation.title
        # ""`) to the end of the file may help.")
        file.writelines('_citation.title  "Elfin"')

def save_pdb(**kwargs):
    '''
    Saves a Bio.PDB.Structure.Structure as a PDB file.

    Args:
    - struct - Bio.PDB.Structure.Structure to be saved.
    - save_path - string file path.
    '''
    struct = kwargs.pop('struct')
    save_path = kwargs.pop('save_path')

    io = Bio.PDB.PDBIO()
    io.set_structure(struct)
    io.save(save_path)

def main():
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()