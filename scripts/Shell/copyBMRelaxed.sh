#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_relaxed"

mkdir -p bm/$outDir
for dir in `ls -d bm/$bmDir`; do 
	for f in `ls bm/$bmDir/$dir/*_0001.pdb`; do 
		cp $f bm/$outDir/$dir"_0001.pdb" 
		break 1			
	done; 
done
