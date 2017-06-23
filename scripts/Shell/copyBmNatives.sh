#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_nat"

mkdir -p bm/$outDir
for dir in `ls -d bm/$bmDir`; do 
	for f in `ls bm/$bmDir/$dir/*.pdb`; do 
		if [[ $f == *"_0001.pdb" ]]; then
			continue
		fi

		cp $f bm/$outDir/$dir".pdb" 
		break 1			
	done; 
done
