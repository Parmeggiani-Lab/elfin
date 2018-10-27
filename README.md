# elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)

Elfin has been split into the following repositories:
 1. [elfin-ui](https://github.com/joy13975/elfin-ui)
 2. [elfin-library](https://github.com/joy13975/elfin-library)
 3. [elfin-data](https://github.com/joy13975/elfin-data) (private)
 4. [elfin-solver](https://github.com/joy13975/elfin-solver)

This repository hosts the v1 branch and data processing code for v2.

Development is currently most active in elfin-ui. Implementation is based on the theories and assumptions stated [here](theories_and_assumptions.md).

# What is elfin

Elfin is a computational protein design tool based on [repeat protein assembly](https://www.sciencedirect.com/science/article/pii/S1047847717301417). 

[Get Started](#2-prerequisites)

Repeat protein assembly means repeat proteins are used as basic building blocks to construct larger proteins that form a 3D structure as close to the user's input description as possible. Credits to Fabio Parmeggiani (UoB), TJ Brunette (UoW), David Baker (UoW), and Simon McIntosh-Smith (UoB).

* UoB: University of Bristol
* UoW: University of Washington

Elfin v2 is an overhaul to add complex design capabilities to the original proof-of-concept ([branch v1](https://github.com/joy13975/elfin/tree/v1)). This resository hosts the data processing scripts, the core loop-closer, and Elfin Front (Blender addon). 

The PDB files of the building blocks are hosted in a [private repository](https://github.com/joy13975/elfin-db). These data have not yet been published by their authors (F. Parmeggiani and Baker lab of the University of Washington), I cannot make them public. If you need access to these data, please contact [Fabio Parmeggiani](https://github.com/parmef) at fabio.parmeggiani@bristol.ac.uk. You can still run Elfin and design proteins without PDB data, but you will not be able to create the full atom model for your design.

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: the handwritten word "Bristol" drawn using protein modules, assembled by Elfin. Visualisation created using [PyMol](https://pymol.org).

### Content:
1. [Project Status](#1-project-status-v2)

2. [Prerequisites](#2-prerequisites)

3. [Setup](#3-setup)

4. [Protein Design UI](#4-protein-design-ui)

5. [Core Solver](#5-core-solver)

6. [Database Preparation](#6-database-preparation)

7. [Stitching Example (v1)](#7-stitching-example-v1)

## 1. Project Status (v2)
- [x] Integrate new hub data
    - [x] Clean, align, extract transformations
- [x] New Maths to more flexibly manipulate complex designs 
- [ ] Reimplement the design synthesis stage
    - [x] Redesign input data format
    - [x] Handle capping
    - [x] Handle multiple networks (1 network = 1 chain)
    - [ ] Handle hub-induced multi-chain network (1 network = N chains)
- [x] UI (Blender addon)
    - [x] Process Pymol .obj files into a Blender library
    - [x] Handle module placement & extrusion
        - [x] Singles
        - [x] Hubs
    - [x] Feasibility validation
        - [x] Collision detection
            - [x] When placing and extruding modules
            - [x After object transform
        - [x] Symmetric hub arm symmetry enforcement
- [x] Handle mixing of module objects with path guides and Key Points
- [x] Handle leeway specification on networks and Key Points
- [ ] Export to stitch.py-readable format
- [x] Create a valid H-shaped design 
- [ ] Elfin Core
- [ ] Solve cyclic multi-chain design
- [ ] Handover
    - [ ] Code documentation
    - [ ] README.md update
        - [ ] Script and directory refactoring
        - [ ] Rewrite installation and running instructions
    - [ ] Simplify installation
        - [ ] Automate python setup
        - [x] Automate installation of Blender addon (install_belfin)
- [ ] Extras
    - [ ] Efficiency optimisation; GPU
    - [ ] Call Elfin Core from Blender - live design?

## 2. Prerequisites
1. [Python 3+](https://www.python.org/downloads/)
2. [VirtualEnv](https://virtualenv.pypa.io/en/stable/)
3. [Blender](https://www.blender.org/)
4. [gcc-5+](https://gcc.gnu.org/)

#### Optional Tools
1. [PyMOL](https://www.pymol.org) for protein visualisation
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for protein relaxation.

## 3. Setup

Run the following:
```
bash <(curl -s https://raw.githubusercontent.com/joy13975/elfin/master/setup_elfin)
```

Note that you will be prompted to enter your Github username and password so as to authenticate for permission to the [elfin-data](https://github.com/joy13975/elfin-data) repo.

## 4. Protein Design UI

See [elfin-ui](https://github.com/joy13975/elfin-ui)

## 5. Core Solver

Example run:
```
./bin/elfin -i ../resources/examples/horns_input.csv -msg 10
```

See [elfin-solver](https://github.com/joy13975/elfin-solver)

## 6. Synthesis Resources

See [elfin-data](https://github.com/joy13975/elfin-data).

## 7. Stitching Example (v1)

If you want to invoke `stitch.py` on a v1 elfin output json file, first convert it to a `stitch.py` readable format. Taking `resources/examples/horns_output.json` as an example. Run at elfin root:
```
. ./activate
./elfinpy/v1_design_convert.py resources/examples/horns_output.json
```

And then you will be able to stitch the v2 output:
```
./elfinpy/stitch.py resources/examples/horns_output.v2.json
```