# elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)

[日本語のご案内はこちら](README_JP.md)

Elfin is suite of tools that helps protein designers build structures by using smaller proteins as building blocks. Elfin does not make any OS/vendor/arch specific assumptions but it is only tested on Linux, MacOS, and WSL.

Elfin is split into the following repositories:
 1. [elfin-ui](https://github.com/joy13975/elfin-ui) - the GUI.
 2. [elfin-library](https://github.com/joy13975/elfin-library) - the public data.
 3. [elfin-data](https://github.com/joy13975/elfin-data) - the private data.
 4. [elfin-solver](https://github.com/joy13975/elfin-solver) - the core.

This main repository hosts data processing scripts for v2. 

There is a v1 branch, which hosts the v1 version of elfin, updated only to be compatible with *some* new scripts. The code there is no longer supported.

Implementation is based on the theories and assumptions stated [here](theories_and_assumptions.md).

# What is elfin?

Elfin is a computational protein design tool suite based on [repeat protein assembly](https://www.sciencedirect.com/science/article/pii/S1047847717301417). 

[Skip To: Get Started](#2-prerequisites).

Repeat protein assembly uses repeat proteins as basic building blocks to construct larger proteins that form a 3D structure as close to the user's input description as possible. Credits to Fabio Parmeggiani (UoB), TJ Brunette (UoW), David Baker (UoW), and Simon McIntosh-Smith (UoB).

* UoB: University of Bristol
* UoW: University of Washington

Elfin v2 is an overhaul to add complex design capabilities to the original proof-of-concept ([branch v1](https://github.com/joy13975/elfin/tree/v1)).

The PDB files of the building blocks are hosted in a [private repository](https://github.com/joy13975/elfin-db). It was requested that these data be kept private before publishing (designed by F. Parmeggiani and Baker lab of the University of Washington). If you need access to these data, please contact [Fabio Parmeggiani](https://github.com/parmef) at fabio.parmeggiani@bristol.ac.uk. You can still run Elfin and design proteins without PDB data, but you will not be able to create the full atom model for your design at the final stage.

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: the handwritten word "Bristol" drawn using protein modules, assembled by Elfin. Visualisation created using [PyMol](https://pymol.org).

### Content
1. [Project Status](#1-project-status)

2. [Prerequisites](#2-prerequisites)

3. [Setup](#3-setup)

4. [Protein Design UI](#4-protein-design-ui)

5. [Autodesign via Machine Learning](#5-autodesign-via-machine-learning)

6. [Protein Data Preprocessing](#6-protein-data-preprocessing)

7. [Creating Output](#7-creating-output)

## 1. Project Status

Functionality is mostly complete for v2, except for some minor left over TODOs noted in issues in respective repositories.

## 2. Prerequisites
#### Required
1. [Python 3+](https://www.python.org/downloads/)
2. [Virtualenv](https://virtualenv.pypa.io/en/stable/)
3. [Blender](https://www.blender.org/)
4. [gcc-5+](https://gcc.gnu.org/)

#### Optional Tools
1. [PyMOL](https://www.pymol.org) for protein visualisation
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for protein relaxation.

## 3. Setup

Run the following:
```Bash
bash <(curl -s https://raw.githubusercontent.com/joy13975/elfin/master/setup_elfin)
```

Note: In order to authenticate for permission to the [elfin-data](https://github.com/joy13975/elfin-data) repo, you will be prompted to enter your Github username and password. If you have not been granted permission, you can just skip this step by hitting enter twice.

## 4. Protein Design UI

See [elfin-ui](https://github.com/joy13975/elfin-ui)

## 5. Autodesign via Machine Learning

See [elfin-solver](https://github.com/joy13975/elfin-solver)

## 6. Protein Data Preprocessing

See [elfin-data](https://github.com/joy13975/elfin-data).

## 7. Creating Output

If you want to invoke `stitch.py` on a v1 elfin output json file, first convert it to a `stitch.py` readable format. Taking `resources/examples/horns_output.json` as an example. Run at elfin root:
```Bash
. ./activate
./elfinpy/v1_design_convert.py resources/examples/horns_output.json
```

And then you will be able to stitch the v2 output:
```Bash
./elfinpy/stitch.py resources/examples/horns_output.v2.json
```

If you want to process solver output entirely in code, below is an example:


```Bash
. ./activate # make sure you're in elfin root and in virtual env
```

Then either in a script or an interactive shell:

```Python
import elfinpy.v1_design_convert as converter
import elfinpy.utilities as utils
import elfinpy.pdb_utilities as pdb_utils
import elfinpy.stitch as stitch

# Normally, you'd need to read a JSON file for the output data from solver.
# Alternatively, create the solver output in Python.

input_json = utils.read_json('resources/examples/horns_output.json')
graphs = converter.v1_to_v2(input_json, 'resources/xdb.json')

# Convert List[graph] to dictionary
graphs_dict = utils.to_dict(graphs)

# Run stitch synth
syn = stitch.Synthesiser(graphs_dict, 
    pdb_dir='./resources/pdb_aligned/',
    cappings_dir='./resources/pdb_relaxed/cappings',
    metadata_dir='./resources/metadata/',
    show_fusion=False,
    disable_capping=False)
struct = syn.run()

as_cif = True

if as_cif:
    pdb_utils.save_cif(struct=struct, path='test.cif')
else:
    # To save as PDB, the chain ID needs to be changed to
    # one that PDB format accepts i.e. a single character.
    # The following line assumes there's only 1 chain. By
    # right in a v1 format there should only be 1 chain 
    # anyway.
    [c for c in struct[0].get_chains()][0].id = 'a'
    pdb_utils.save_pdb(struct=struct, path='test.pdb')
```