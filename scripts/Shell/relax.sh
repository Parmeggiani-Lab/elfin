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

defaultMaxCycles=200
defaultLocal="yes" 					# or no
defaultVariant="omp" 				# or mpi
defaultRelease="linuxgccrelease"	# or macosclangrelease 
defaultWrapper="" 					# or mpirun

maxCycles=${maxCycles:-$defaultMaxCycles}
local=${local:-$defaultLocal}
variant=${variant:-$defaultVariant}
release=${release:-$defaultRelease}
wrapper=${wrapper:-$defaultWrapper}

if [[ "$variant" == "mpi" ]]; then
	wrapper=$wrapper" mpirun"
fi

cmd="$wrapper relax.$variant.$release -overwrite -s $input -out:path:score $outDir -out:file:scorefile $scOutput -out:path:pdb $outDir -default_max_cycles $maxCycles"

if [[ "$local" == "yes" ]]; then
	$cmd
else
	sbatch -A other -p cpu -N 1 --ntasks-per-node=16 --wrap="$cmd"
fi
