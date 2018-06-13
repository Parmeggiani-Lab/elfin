#/bin/bash

if [[ $# -lt 1 ]]; then
	echo 'Usage: minAll <dir>'
	exit
fi

bmDir="$1"

for f in `ls *.pdb`
do
	pdbOutput="$bmDir/minimised/${f/.pdb/_0001.pdb}"

	if [ -f $pdbOutput ]; then
		echo "$f already minimised - skipping"
	else
		input="$bmDir/$f"
		scOutput="$bmDir/${f/pdb/min}.sc"

    	echo "Minimising $f"
	    minimize.default.macosclangrelease -overwrite \
	        -s $input \
	        -out:path:score scores -out:file:scorefile $scOutput \
	        -out:path:pdb minimised
	fi
done
