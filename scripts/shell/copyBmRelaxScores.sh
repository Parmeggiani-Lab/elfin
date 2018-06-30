#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_relax_sc"

mkdir -p bm/$outDir
for dir in `ls -d bm/$bmDir`; do 
	for f in `ls bm/$bmDir/$dir/*_comp.sc`; do 
		cp $f bm/$outDir/$dir"_comp.sc" 
		break 1			
	done; 
done
