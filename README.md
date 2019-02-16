# elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)

[日本語のご案内はこちら](README_JP.md)

Elfin is suite of tools that helps protein designers build structures by using smaller proteins as building blocks. Elfin does not make any OS/vendor/arch specific assumptions but it has only been tested on Linux, MacOS, and WSL.

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

[Skip To: Usage](#usage).

Repeat protein assembly uses repeat proteins as basic building blocks to construct larger proteins that form a 3D structure as close to the user's input description as possible. Credits to Fabio Parmeggiani (UoB), TJ Brunette (UoW), David Baker (UoW), and Simon McIntosh-Smith (UoB).

* UoB: University of Bristol
* UoW: University of Washington

Elfin v2 is an overhaul to add complex design capabilities to the original proof-of-concept ([branch v1](https://github.com/joy13975/elfin/tree/v1)).

The PDB files of the building blocks are hosted in a [private repository](https://github.com/joy13975/elfin-db). It was requested that these data be kept private before publishing (designed by F. Parmeggiani and Baker lab of the University of Washington). If you need access to these data, please contact [Fabio Parmeggiani](https://github.com/parmef) at fabio.parmeggiani@bristol.ac.uk. You can still run Elfin and design proteins without PDB data, but you will not be able to create the full atom model for your design at the final stage.

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: the handwritten word "Bristol" drawn using protein modules, assembled by Elfin. Visualisation created using [PyMol](https://pymol.org).

## Project Status

Functionality is mostly complete for v2, except for some minor leftover TODOs noted as issues in respective repositories.

## Usage
### Content
   1. [Prerequisites](#1-prerequisites)
   2. [Setup](#2-setup)
   3. [Workflow](#3-design-workflow)
   4. [Creating Output for v1](#4-creating-output-for-v1)
   5. [Scripting](#5-scripting)

### 1. Prerequisites
Firstly install the following software:
#### Required
1. [Python 3+](https://www.python.org/downloads/) for data processing scripts
2. [Virtualenv](https://virtualenv.pypa.io/en/stable/) ^

#### Optional Tools
1. [PyMOL](https://www.pymol.org) for protein visualisation
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for protein optimization

### 2. Setup

Run the following command which calls the auto setup script:
```Bash
bash <(curl -s https://raw.githubusercontent.com/joy13975/elfin/master/setup_elfin)
```

Note: In order to authenticate for permission to the [elfin-data](https://github.com/joy13975/elfin-data) repo, the script will ask you to enter Github username and password. If you have not been granted permission, you can skip this step by hitting enter twice.

### 3. Design Workflow

#### Create Geometry
The workflow begins with drawing out the shape you would like to build using proteins. This is done via [elfin-ui](https://github.com/joy13975/elfin-ui), where corresponding documentation is  available. 

After drawing the specification in Blender using elfin's plugin, export it (elfin-ui command: #exp) to a JSON file.

#### Autodesign
Next, use [elfin-solver](https://github.com/joy13975/elfin-solver) to auto-design the target geometry. The solver outputs another JSON file with design solutions.

#### Fixup
In Blender, open a blank file and import (elfin-ui command: #imp) the solver's output.

Fix up the solution if needed be (perhaps by closing any gaps or joining separated networks).

Export (#exp) the solution again but this time, ensure no path guide objects are present because the next stage will not accept a JSON with path guide objects.

#### Export as PDB/CIF
Lastly, run the following command:

```
. ./activate  # Activates venv
stitch.py <PATH_TO_YOUR_EXPORTED_SOLUTION_JSON>
```

Note that ```. ./activate``` only needs to be run when `venv` is not active.

#### For those who wish to do/redo data preprocessing:

Protein data has already been preprocessed and hosted in elfin-data, so for most people this step is not needed. If new data has been added to the module database or the preprocessing method has changed, then you may wish the redo the data preprocessing. 

See [elfin-data](https://github.com/joy13975/elfin-data)(private).

### 4. Creating Output for v1

This is no longer supported due to a breaking change in the `stitch.py`. There should be no need to do this anymore since elfin-solver v2 supports the same functionality for v1.

### 5. Scripting

You can use elfin's data processing classes in your own script or interactive shell. Below is an example:

```Python
import elfinpy.utilities as utils
import elfinpy.pdb_utilities as pdb_utils
import elfinpy.stitch as stitch

struct = stitch.Stitcher(
    spec=utils.read_json('resources/examples/half_snake_2x1h_deposit_test.json'),
    xdb=utils.read_json('resources/xdb.json'),
    pdb_dir='./resources/pdb_aligned/',
    cappings_dir='./resources/pdb_relaxed/cappings',
    metadata_dir='./resources/metadata/',
    show_fusion=False,
    disable_capping=False,
    skip_unused=False).run()

as_cif = True

if as_cif:
    pdb_utils.save_cif(struct=struct, path='test.cif')
else:
    cid = 'A'
    for c in struct[0].get_chains():
        if cid > 'Z':
            raise ValueError('Too many chains for PDB format')
            # You might want to extend number of chains by using lowercase
            # alphabets.

        c.id = cid
        cid = chr(ord(cid) + 1)

    pdb_utils.save_pdb(struct=struct, path='test.pdb')
```