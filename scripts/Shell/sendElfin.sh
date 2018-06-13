	#!/bin/bash

: ${jobName=elfin};
: ${queueName="batch"};
: ${ppn="256:knl"};
: ${runScript="./runElfin.sh"};
: ${walltime="00:30:00"};
: ${execFile='./src/cpp/bin/elfin'};
: ${configFile='elfinConfig.json'};
: ${outputDir="output"};
: ${extraArgs=''};

cmd="qsub -N $jobName \
	-q $queueName \
	-joe -o $outputDir/log \
	-l nodes=1:ppn=$ppn,walltime=$walltime \
	-v execFile=$execFile,configFile=$configFile,outputDir=$outputDir,extraArgs=$extraArgs \
	 $runScript"

mkdir -p $outputDir
echo $cmd > $outputDir/cmd
$cmd
