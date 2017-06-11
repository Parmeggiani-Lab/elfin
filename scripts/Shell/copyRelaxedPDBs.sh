#!/bin/bash

inputDir=${inputDir:-"./res/preprocessed/"}
outputDir=${outputDir:-"./res/relaxed/"}

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