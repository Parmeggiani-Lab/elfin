#/bin/bash

# Or mpi, if available
variant=default

# For Macs: macosclangrelease
releaseString=linuxgccrelease

for f in `ls output_*/*_Synth_*.pdb`
do
	pdbOutput="${f/\.pdb/_0001.pdb}"

	if [ -f $pdbOutput ]; then
		echo "$f already minimised - skipping"
	else
		scOutput="${f/\.pdb/_min.sc}"
		outDir=`dirname $f`
    		echo "Minimising $f; Output: $pdbOutput; MinScore: $scOutput"

		#exit
		minimize.$variant.$releaseString -overwrite \
	        	-s $f \
	        	-out:path:score $outDir -out:file:scorefile $scOutput \
		        -out:path:pdb $outDir
	fi
done
