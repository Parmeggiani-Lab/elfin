#!/bin/bash
trap "exit" INT

bmDir=${bmDir:-"../../bm/l20/"}
gps=${gps:-262144}

export OMP_SCHEDULE=dynamic,16

mkdir -p done
for f in `ls $bmDir/*.json`
do
	mkdir -p output

	cmd="./bin/elfin -i $f -gps $gps"
	echo Running cmd: $cmd
	$cmd

	bmName=${f/\.json//}
	leading=`dirname $f`'//'
	bmName=${bmName/$leading//}
	dest=done/$bmName

	echo Moving output to $dest
	mv output done/$bmName
done