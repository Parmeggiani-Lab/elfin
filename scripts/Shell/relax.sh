#!/bin/bash

trap "exit" INT

if [ $# -lt 1 ]; then
	echo "Usage: relax.sh <input_pdb> <max_cycles=200>"
	exit
fi

input="$1"
maxCycles="$2"

outDir=`dirname $input`
scOutput="${input/\.pdb/_relax.sc}"

maxCycles=${maxCycles:-200}
local=${local:-"no"}
variant=${variant:-"mpi"}
release=${release:-"linuxgccrelease"}
wrapper=${wrapper:-""}

if [[ "$variant" == "mpi" ]]; then
	wrapper=$wrapper" mpirun"
fi

cmd="$wrapper relax.$variant.$release -overwrite -s $input -out:path:score $outDir -out:file:scorefile $scOutput -out:path:pdb $outDir -default_max_cycles $maxCycles"

if [[ "$local" == "yes" ]]; then
	$cmd
else
	sbatch -A other -p cpu -N 1 --ntasks-per-node=16 --wrap="$cmd"
fi
