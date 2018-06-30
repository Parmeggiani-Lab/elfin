#!/bin/bash

trap "exit" INT

if [ $# -lt 1 ]; then
	echo "Usage: relax.sh <input_pdb> <max_cycles=200> <overwrite=no>"
	exit
fi

input="$1"
maxCycles="$2"
overwrite="$3"

outDir=`dirname $input`
scOutput="${input}_relax.sc"

defaultMaxCycles=200
defaultLocal="yes" 					# or no
defaultVariant="omp" 				# or mpi
defaultRelease="linuxgccrelease"	# or macosclangrelease 
defaultWrapper="" 					# or mpirun

maxCycles=${maxCycles:-$defaultMaxCycles}
local=${local:-$defaultLocal}
overwrite=${overwrite:-"no"}
variant=${variant:-$defaultVariant}
release=${release:-$defaultRelease}
wrapper=${wrapper:-$defaultWrapper}

if [[ "$variant" == "mpi" ]]; then
	wrapper=$wrapper" mpirun"
fi

cmd="$wrapper relax.$variant.$release -overwrite -s $input -out:path:score $outDir -out:file:scorefile $scOutput -out:path:pdb $outDir -default_max_cycles $maxCycles"
echo "cmd="$cmd
if [ ! -f ${input/.pdb/_0001.pdb} ] || [[ "$overwrite" == "yes" ]]; then
	if [[ "$local" == "yes" ]]; then
		$cmd
	else
		sbatch -A other -p cpu -N 1 --ntasks-per-node=16 --wrap="$cmd"
	fi
fi
