#!/bin/bash

trap "exit" INT

config=${config:-"src/GA/config.json"}
xdb=${xdb:-"res/xDB.json"}
gps=${gps:-524288}
local=${local:-"no"}
dir=${dir:-"bm/l10"}
outdir=${outdir:-$dir}

export OMP_SCHEDULE=dynamic,16

for f in `ls $dir/*.json`
do	
	filename=`basename $f`
	myOutdir=$outdir'/'${filename/\.json/}
	mkdir -p $myOutdir

	cmd="./src/GA/bin/elfin -i $f -c $config -gps $gps -x $xdb -o $myOutdir"
	echo CMD is $cmd
  	if [[ "$local" == "yes" ]]; then
#		echo local
		$cmd
	else
#		echo remote
        	sbatch -A other -p cpu -N 1 --ntasks-per-node=16 --wrap="$cmd"
	fi
done

