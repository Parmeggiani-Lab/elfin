#!/bin/bash

inputDir=${1:-"./res/preprocessed_modules/"}
outputDir=${2:-"./res/aligned_relaxed_modules/"}

mkdir -p $outputDir'/single/'
mkdir -p $outputDir'/pair/'

for f in `ls $inputDir/single/*_0001.pdb`; do
	bn=`basename $f`
	cp $f $outputDir'/single/'${bn/_0001/}
done

for f in `ls $inputDir/pair/*_0001.pdb`; do
	bn=`basename $f`
	cp $f $outputDir'/pair/'${bn/_0001/}
done