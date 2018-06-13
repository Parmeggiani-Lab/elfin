#!/bin/bash

default_input_pdbs="./res/preprocessed_modules/*/*.pdb"
input_pdbs=${1:-$default_input_pdbs}
env_vars=${2:-""} 

if [[ "$1" == "-h" ]]; then
	echo "Usage: dbRelaxList.sh <input_pdbs=\""${default_input_pdbs}"\"> <env_vars=none>"
	echo "	Make sure to quote the input_pdbs argument!"
	exit
fi

echo 'trap "exit" INT'

ls ${input_pdbs} | xargs -I{} echo ${envs}' ./scripts/Shell/relax.sh {}'