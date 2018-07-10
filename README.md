# Elfin v2 [![Build Status](https://travis-ci.com/joy13975/elfin.svg?branch=master)](https://travis-ci.com/joy13975/elfin)
Elfin is a computational protein design tool based on [repeat protein assembly](https://www.sciencedirect.com/science/article/pii/S1047847717301417). 

To jump straight in, skip to [Prerequisites](#1-prerequisites).

Repeat protein assembly means repeat proteins are used as basic building blocks to construct larger proteins that form a 3D structure as close to the user's input description as possible. Credits to Fabio Parmeggiani (UoB), TJ Brunette (UoW), David Baker (UoW), and Simon McIntosh-Smith (UoB).

* UoB: University of Bristol
* UoW: University of Washington

Elfin v2 is an overhaul to add complex design capabilities to the original proof-of-concept ([branch v1](https://github.com/joy13975/elfin/tree/v1)). This resository hosts the data processing scripts, the core loop-closer, and Elfin Front (Blender addon). 

The PDB files of the building blocks are hosted in a [private repository](https://github.com/joy13975/elfin-db). These data have not yet been published by their authors (F. Parmeggiani and Baker lab of the University of Washington), I cannot make them public. If you need access to these data, please contact [Fabio Parmeggiani](https://github.com/parmef) at fabio.parmeggiani@bristol.ac.uk. You can still run Elfin and design proteins without PDB data, but you will not be able to create the full atom model for your design.

![alt tag](resources/diagrams/ProteinBristol.png)
Figure 1: the handwritten word "Bristol" drawn using protein modules, assembled by Elfin. Visualisation created using [PyMol](https://pymol.org).

### Content:
1. [Project Status](#1-project-status)

2. [Prerequisites](#2-prerequisites)

3. [Optional Tools](#3-optional-tools)

4. [Repository Setup](#4-repository-setup)

5. [Python Setup](#5-python-setup)

6. [Core Setup](#6-core-setup)

7. [Usage](#7-usage)

8. [Protein Design UI](#8-protein-design-ui)

9. [Database Preparation](#9-database-preparation)

10. [Design Verification](#10-design-verification)

## 1. Project Status

### v1
- [x] Allows user to draw simple shapes in Matlab to describe the protein they want to make
- [x] Accepts input as an array of 3D points specified by either a JSON or CSV file
- [x] Scripts for generating database (in case of new modules)
- [x] Produces quality protein designs
- [x] Several tunable parameters in the GA (Core)
- [x] Scripts for design verification (against relaxed versions) with the help of Rosetta 

### v2
- [x] Integrate new hub data
    - [x] Clean, align, extract transformations
- [x] New Maths to more flexibly manipulate complex designs 
- [.] Reimplement the design synthesis stage
    - [x] Redesign input data format
    - [x] Handle capping
    - [x] Handle multiple networks (1 network = 1 chain)
    - [.] Handle hub-induced multi-chain network (1 network = N chains)
- [.] UI (Blender addon)
    - [x] Process Pymol .obj files into a Blender library
    - [.] Handle module placement & extrusion
        - [x] Singles
        - [.] Hubs
- [ ] Handle mixing of module objects with path guides and Key Points
- [ ] Handle leeway specification on networks and Key Points
- [ ] Blend & Export to stitch.py-readable format
- [ ] Create a valid H-shaped design 
- [ ] Elfin Core
- [ ] Solve cyclic multi-chain design
- [.] Handover
    - [.] Code documentation
    - [.] README.md update
        - [.] Script and directory refactoring
        - [ ] Rewrite installation and running instructions
    - [.] Simplify installation
        - [ ] Automate python setup
        - [x] Automate installation of Blender addon (install_belfin)
- [ ] Extras
    - [ ] Efficiency optimisation; use GPU if it helps
    - [ ] Call Elfin Core from Blender - live design?

## 2. Prerequisites
1. [Python 2.7+](https://www.python.org/downloads/release/python-279/)
2. [VirtualEnv](https://virtualenv.pypa.io/en/stable/)
3. [Blender](https://www.blender.org/)
4. gcc-5+
5. 7z, either its GUI or CLI will do

## 3. Optional Tools
1. [PyMOL](https://www.pymol.org) for protein visualisation
2. [Rosetta](https://www.rosettacommons.org/software/license-and-download) for protein relaxation.

## 4. Repository Setup

```
git clone https://github.com/joy13975/elfin.git
cd elfin                                                
```

## 5. Python Setup

```
# cd to repo root
virtualenv .venv
. ./activate
pip install -r requirements.txt
```

## 6. Core Setup

First you need GCC/G++ 5 and above. On Macs for instance, you can get OpenMP-enabled GCC-6 by simply installing from Homebrew:

```
brew install gcc --without-multilib
```

The installation could take a while. After this is done, you should get commands ```gcc-6``` and ```g++-6```. You can choose to make these your default compilers by overwriting the symbolic links at ```/usr/local/bin/gcc``` and ```/usr/local/bin/g++```:
```
ln -Fs `which gcc-6` /usr/local/bin/gcc
ln -Fs `which g++-6` /usr/local/bin/g++
```

For troubleshooting on MacOS, please refer to [this Stackoverflow post](https://stackoverflow.com/questions/35134681/installing-openmp-on-mac-os-x-10-11). 

If your compiler complains about ```omp_get_initial_device``` not declared, that's because your OpenMP version is too old. Check with ``` echo |cpp -fopenmp -dM |grep -i open```; only versions above 201511 seem to define this function. This is not vital to the application so if you really do not want to install a newer compiler then you may comment out that erring line.

Then install the jutil submodule and get ready to compile.

```
git submodule init
git submodule update --init --force --remote            
cd celfin                                                   
make
```

Notes:
 - Specify the compiler of your choice by e.g. for clang++: ```make CXX=clang++```.
 - For clang, you will need to specify the C++ standard library include path, which depends on where you installed GCC. See ./celfin/Makefile for details (the INCS variable). You will also need to include libiomp and library load paths (again see Makefile).
 - To speed up the compilation, add the ```-j``` flag with the number of cores on your system.
 - To build without OpenMP, you can specify ```make OMP=no```

## 7. Usage
Once you have compiled Elfin Core successfully, you can test run it with:
```
# cd to celfin
./bin/elfin                                             
```

The default configuration is at ```./celfin/config.json``` and uses the ```./bm/l10/1696.json``` example input. This demonstrates that the Core GA finds the zero score solution i.e. the original module constituents of the input. Execution outputs will be stored in ```./celfin/output/``` in JSON format (files named after their memory addresses). You can specify the output directory using the ```-o``` flag. Even if you force Elfin to stop using ctrl-c, the most recent best candidates will be saved.

To test a different example, you may try:

```
./bin/elfin -i ../bm/l20/mqj2.json
```

All examples in ```./bm/l10```, ```./bm/l20``` and ```./bm/l30``` should reach score zero before the maximum number of iterations (default=1000) is exhausted. However, because the random number generator might behave differently on your machine, in some rare cases Elfin might not reach score zero. If this happens, try using larger population sizes e.g. using the ```-gps``` flag.

To get help:
```
./bin/elfin -h
```

A typical Elfin run with specific population size and input file would be:
```
./bin/elfin -gps <POPULATION_SIZE> -i <INPUT_FILE>
```

Notes:
 - Output JSONs are named after their sorted memory addresses. The lower value the hex name has, the better solution it is (lower score).
 - Command-line arguments will override arguments specified in the configuration file.

## 8. Protein Design UI

First, decompress the Blender protein model library:
```
# cd to resources
7z x elfin_front_library.7z
```

Next, make sure Blender is installed. Then invoke the addon installation script:
```
# cd to repo root and make sure venv is activated
install_belfin
```

If you're on MacOS, the script will probably guess your Blender directory and let you choose a version. Always choose the latest version.

If you're on Linux or WSL, the script will ask you to provide your Blender directory. Again, choose the latest version.

After the script is done symlinking, open Blender and go to ```File->User Preferences->Add-ons```, and search for "Elfin". Tick the item named "Elfin: Elfin Front". Hit "Save User Settings" so that this gets remembered for future launches.

If the above went successfully, your left-hand-side vertical bar should have a new "Elfin" panel. If the bar is not visible, toggle it with <kbd>T</kbd>. For now, the buttons in this panel are for debugging, so a non-developer user should not find any business here.

Now you're ready to test Elfin Front's features.

Inside the 3D view port, make sure that you're in Object mode and viewing Solids. When you hit <kbd>Enter</kbd> you will get a prompt for typing in a command. Elfin Front currently implements the following commands:

 - Place a module
 - Extrude N (add a module to the nterm)
 - Extrude C (add a module to the cterm)
 - (Re)load module libary
 - (Re)load xdb

When selecting a protein module from the given list, you can use ```.``` (the period mark) to indicate the beginning of the module name and the end. This is not compulsory but it helps filter your module of choice.

To delete a placed module, hit <kbd>Del</kbd> or the equivalent on MacOS. Please be warned that a lot of Elfin Front's functionalities are still being developed and tested, so nothing is guaranteed.

Currently the only way to save your design is through the default Save, as a .blend file.

To be updated with more implementation progress.

## 9. Database Preparation

**Obtain and decompress elfin-db.zip**
```
# cd to repo root
. ./activate
cd resources
git clone https://github.com/joy13975/elfin-db.git 
```

Read [here](#access-to-elfin-db) about elfin-db access permission.

Obtain 7z. For Ubuntu:
```
sudo apt-get install p7zip-full
```
And then decompress:

```
# cd to elfin-db inside resources
7z x elfin-db/db.7z
```


**Preprocess and relax module PDBs**
```
# cd to repo root
preprocess.py
create_relax_list > relax_list.temp

```

Now, if you're on a Linux cluster **that has MPI-enabled Rosetta installed** and you have the ```sbatch``` command for your job queuing system, you can simply invoke:

```
sh relax_list.temp
```

If you are on your own computer or some other different environment, you need to modify the environment variables that specify what flavour of Rosetta and queuing system you use. For instance, on my Macbook Pro (without MPI Rosetta and no queuing system) I use:
```
local=yes variant= release=macosclangrelease sh relax_list.temp
```

 - The ```local``` variable tells the ```./scripts/Shell/relax.sh``` script whether or not to use ```sbatch```.
 - The ```variant``` variable tells which build variant (e.g. ```mpi``` or ```none``` if default) of Rosetta to use. 
 - The ```release``` variable tells which release version to use. This depends on your OS, and which build version you want to use (e.g. could be the debug version).

This relaxation can take quite a while simply because the process is computationally intensive. It is therefore strongly recommended that you do this on a compute cluster that lets you spread workload across many machines. 

**Generating xdb**

After relaxation, you should see a new PDB suffixed with ```_0001.pdb``` for each PDB in ```./resources/pdb_prepped/```. These are the relaxed versions of the original preprocessed PDBs. Now copy these into a separate folder for the next step:
```
cp_relaxed_pdbs
```

Lastly, generate the xdb by invoking: 
```
dbgen.py
```

This script aligns each pair's reference frames to the reference frame of their first single module (e.g. ```D4-D4_j1_D64``` would be superimposed onto ```D4```), producing the ```./resources/pdb_aligned/``` folder. This folder contains aligned singles and pairs that will be used to reconstruct Elfin's output structures in the final stage of design. Geometric relationships will also be calculated to produce a ```./resources/xdb.json``` "pair relationship database".

## 10. Design Verification
After running a shape design through Elfin Core, you should find output solutions in ```./celfin/output/```. Recall that smaller hexadecimal values in the file name corresponds to better solutions (they're just memory address of a sorted array). Here we shall use the default input as example.

```
for f in `ls output`; do echo $f; break; done #print best solution file name
# take for instance 0x10e63c000.json
# go back to repo root       
```
Inspect the output file to see what it contains (score and module/node names). The score is the Kabsch "simularity" score between the solution and the input shape. The lower this score is, the better the solution.

**PDB Reconstruction**

To be re-written.

**Visualisation**

At this point, you can already visualise the design outputs in PyMOL by loading the synthesised PDB. You can also plot the solution centre-of-mass points (stored in the CSV output mentioned above). In PyMOL: load the plotting script from ```./scripts/PyMolUtils/LineUtils.py```. Then, type the following in the command window:
```
draw_csv('path/to/your/csv/file')                       #to plot any CSV 3D point file
draw_json('path/to/your/json/file')                     #to plot JSON files that have a 'coms' 3D point array
```

More plotting options can be found by reading the script ```./scripts/PyMolUtils/LineUtils.py```. Note that the default plotting color is black, so you must do ```bg white``` to see them properly. Also note that input specification CSV files can be plotted using the same command.

**Calculating Global RMSD**

To compute the RMSD, we first need to relax the result PDB. This is the same as when we relaxed databse PDBs:

```
relax </path/to/PDB/for/relaxation>
```

Apply the same environment variable changes as needed (mentioned in section 6).
The relaxation of larger structures will take longer (much much longer than DB modules). Using our previous example (on ```./bm/l10/1696.json```) we should yield a file called ```./celfin/output/0x10e63c000_Synth_Main_s1.0_l10_0001.pdb``` and ```./celfin/output/0x10e63c000_Synth_Main_s1.0_l10_0001_relax.sc```. 

The first file is the relaxed PDB structure, and the second is the score file that contains information about how different the structures are before and after relaxation. To get a Rosetta RMSD estimate (global), open the score file and look at the second column. In this case, I got 2.424.

In the score file the first row is the column headers. Each subsequent row denotes values obtained from each subsequent Rosetta transformation of the structure in question. It is useful to divide the global RMSD by the number of modules in the design to get a more unbiased perspective of RMSD per module. This is because lever effect is obviously much more pronounced in longer/larger protein designs.

**Calculating Windowed RMSD**

Ensure that you have a relaxed result PDB first. To compute the windowed RMSD, invoke:
```
rmsd.py <solutionDir> <minimisedDir> <windowLen=300> <overlapRatio=0.50> <warnThreshold=5.0>
```

Before you can run this script, copy the relaxed (or minimised) PDBs you want to compare into a different folder - one which the script shall treat as ```<minimisedDir>```. For each relaxed/minimised PDB in this folder, the script will compute a windowed RMSD with respect to their original structure (should exist in ```<solutionDir>```) before the Rosetta transformation.

For instance if we run this on our example output:
```
mkdir ./tmp                                             #make a temporary folder in repo root
cp ./celfin/output/0x10e63c000_Synth_Main_s1.0_l10_0001.pdb ./tmp 
./scripts/Python/RMSDStat.py ./celfin/output ./tmp          #invoke the windowed RMSD script
```

The output indicated:
```
avg: 1.37087369713 min 0.979476939189 max 1.65581158797
```

Which is consistent with the result in my thesis (for benchmark ```./bm/l10/1696.json```).