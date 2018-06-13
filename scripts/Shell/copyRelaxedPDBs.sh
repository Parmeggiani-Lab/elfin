#!/bin/bash

inputDir=${1:-"./resources/pdb_preppd/"}
outputDir=${2:-"./resources/pdb_relaxed/"}

mkdir -p $outputDir'/singles/'
mkdir -p $outputDir'/doubles/'

for f in `ls $inputDir/singles/*_0001.pdb`; do
	bn=`basename $f`
	cp $f $outputDir'/singles/'${bn/_0001/}
done

for f in `ls $inputDir/doubles/*_0001.pdb`; do
	bn=`basename $f`
	cp $f $outputDir'/doubles/'${bn/_0001/}
done