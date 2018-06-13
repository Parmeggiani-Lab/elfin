#/bin/bash

# Or mpi, if available
variant=default

# For Macs: macosclangrelease
releaseString=linuxgccrelease

for f in `ls output_*/*_Synth_*.pdb`
do
	scOutput="${f/\.pdb/_comp.sc}"

	if [ -f $scOutput ]; then
		echo "$f already scored - skipping"
	else
		minPdb="${f/\.pdb/_0001.pdb}"

		echo "Scoring $f; MinPDB: $minPdb; Output: $scOutput"

		exit
		score.$variant.$releaseString \
	        	-s $minPdb \
	        	-in:file:native $f \
	        	-out:file:scorefile $scOutput
	fi
done
