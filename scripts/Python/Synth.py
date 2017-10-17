#!/usr/bin/env python

import glob
import argparse
from utils import *

def main():
    ap = argparse.ArgumentParser(description='Generate PDB from node information');
    ap.add_argument('specFile')
    ap.add_argument('--outFile', default='')
    ap.add_argument('--movieMode', dest='movieMode', action='store_true')
    ap.add_argument('--no-movieMode', dest='movieMode', action='store_false')
    ap.add_argument('--singlesDir', default='res/aligned/single/')
    ap.add_argument('--pairsDir', default='res/aligned/pair/')
    ap.add_argument('--xdbPath', default='res/xDB.json')
    ap.set_defaults(movieMode=False)

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)
        
    args = ap.parse_args()

    specExt = args.specFile[args.specFile.rfind('.'):]

    if specExt == '.json':
        spec = readJSON(args.specFile)
        targetLen = len(spec['nodes'])
    else:
        print 'Unknown spec file type: {}'.format(specExt)
        exit()

    scale = 1.0
    if args.outFile == '':
        args.outFile = args.specFile

    dotIndex = args.outFile.rfind('.')
    if dotIndex == -1:
    	outExt = ''
    	args.outFile = args.outFile + suffixPdb(
		            'Synth', 
		            'Main',
		            scale, 
		            targetLen)
    else:
		outExt = args.outFile[dotIndex:]
		args.outFile = args.outFile.replace(outExt, suffixPdb(
		        'Synth', 
		        'Main',
		        scale, 
		        targetLen))

    xDB = readJSON(args.xdbPath)
    _, comShape = makePdbFromNodes(
        xDB, 
        spec['nodes'], 
        args.pairsDir,
        args.singlesDir,
        args.outFile, 
        movieMode=args.movieMode
    )
    
    csvOutFile = args.outFile.replace('.pdb', '.csv')
    saveCSV(comShape, csvOutFile)

if __name__ == '__main__':
    main()