#!/bin/bash

: ${gsVer=1}
: ${startGridId=0} #includsive
: ${endGridId=0} #inclusive
: ${queueName=teaching}
: ${ppn=16}
: ${gsConfDir=./gsConfigsV$gsVer/}
: ${gsOutDir=./gsOutV$gsVer/}

prog="./src/GA/bin/elfin"
jobCount=0
for gridId in $(seq $startGridId $endGridId); do
	#Each config ID can have multiple benchmarks attached
	bmConfigs=($(ls $gsConfDir/*_"$gridId"_*\.json))

	for bmConfig in "${bmConfigs[@]}"; do
		cmd="$prog -s $bmConfig"
	        bmOutDir=${bmConfig/$gsConfDir/$gsOutDir}
        	bmOutDir=${bmOutDir/\.json/'/'}
        	cmdDir=$bmOutDir/cmd

	        mkdir -p $bmOutDir

	        echo "cd $HOME/src/elfin" > $cmdDir
	        echo $cmd >> $cmdDir

#	        echo cmd is $cmd
#	        echo bmOutDir is $bmOutDir

		jobCount=$((jobCount + 1))
		echo Sending job "#$jobCount ("$bmConfig")"
        	qsubCmd="qsub -N elfin_omp -joe -o $bmOutDir/log -q $queueName -l nodes=1:ppn=$ppn,walltime=00:30:00 ""$cmdDir"
#		echo qsubCmd: $qsubCmd
		$qsubCmd
	done
done
