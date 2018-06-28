#!/bin/bash

if [[ "$#" -lt 2 ]]; then
	echo 'copyRelaxedPDBs.sh <inputDir> <outputDir>'
	exit
fi

inputDir=$1
outputDir=$2

mkdir -p $outputDir

for f in `ls $inputDir/*_0001.pdb`; do
	bn=`basename $f`
	cp $f $outputDir'/'${bn/_0001/}
done
