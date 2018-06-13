#/bin/bash

if [[ $# -lt 1 ]]; then
	echo 'Usage: scoreAll <dir>'
	exit
fi

bmDir="$1"

for f in `ls *.pdb`
do
	scOutput="$bmDir/scores/${f/pdb/comp}.sc"

	if [ -f $scOutput ]; then
		echo "$f already scored - skipping"
	else
	    pdbOutput="$bmDir/minimised/${f/\.pdb/_0001.pdb}"
	    native="$bmDir/$f"

	    echo "Scoring $native"
	    score.default.macosclangrelease \
	        -s $input \
	        -in:file:native $native \
	        -out:file:scorefile $scOutput
	fi
done
