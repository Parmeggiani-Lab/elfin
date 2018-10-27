#!/usr/bin/env python2

#
# This is a PyMol extension script to load all PyMol extensions in the same
# directory.
#

def main():
    """main"""
    raise RuntimeError('This module should not be executed as a script')

if __name__ =='__main__': 
    main()

in_pymol = False
try:
    import pymol
    in_pymol = True
except ImportError as ie:
    main()

if in_pymol:
    from pymol import cmd
    import sys, os, glob

    curr_dir = os.getcwd()
    sys.path.append(os.path.join(curr_dir, os.pardir, os.pardir)) # for elfinpy

    import importlib
    ext_exclusion = ['__init__.py', 'load_all_extensions.py']
    extensions = [
        py for py in glob.glob(os.getcwd() + '/*.py') \
        if os.path.basename(py) not in ext_exclusion \
    ]

    for ext in extensions:
        cmd.load(ext)

    print('All Extensions Loaded')