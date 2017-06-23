#!/bin/bash

bmDir=${bmDir:-"l10"}

for file in `ls $PWD/bm/$bmDir"_nat"/*.pdb`; do
	echo reinit

	echo load $file
	fileDir=`dirname $file`
	fileBase=`basename $file`
	relF=$PWD/bm/$bmDir"_relaxed"/${fileBase/\.pdb/_0001.pdb}
	echo load $relF

	echo disable
	echo hide everything

	obj=${fileBase/\.pdb/}
	relObj=$obj"_0001"

	echo load $PWD/src/PyMolUtils/LineUtils.py
	echo load $PWD/src/PyMolUtils/Compare.py

	specExt='csv'
	specFile=$PWD/bm/$bmDir/$obj.$specExt
	if [ -f "$specFile" ]; then
		echo "draw_"$specExt"('"$specFile"', 1, 2, True)"
	else
		specExt='json'
		specFile=$PWD/bm/$bmDir/$obj.$specExt
		echo "draw_"$specExt"('"$specFile"', 1, 2, True)"
	fi

	solCsv=`ls ${specFile/\.$specExt/}/*.csv`
	rotStr=`./src/Python/GenPymolSpecSolRot.py $specFile $solCsv`
    echo "cmd.transform_selection('"$obj"', $rotStr, homogenous=0)"

	echo bg white
  	echo set depth_cue, 1

	# for chars
	# echo "set_view (\
	#      0.792623580,    0.107972473,   -0.600072622,\
	#      0.000102258,    0.984166384,    0.177220225,\
	#      0.609710038,   -0.140525863,    0.780067921,\
	#      0.000217728,   -0.001177430, -1687.383666992,\
	#     -6.746949673,   93.794700623,  -40.651931763,\
	#   	1256.293457031, 2118.493164062,  -20.000000000 )"
	# For weird shapes
	echo "set_view (\
		 0.682993114,   -0.467070848,    0.561571896,\
		 0.001815517,    0.769914448,    0.638143957,\
		-0.730422139,   -0.434827507,    0.526694298,\
		 0.000000000,    0.000000000, -1833.382080078,\
		 4.250000000,    6.750000000,    0.000000000,\
		1402.282592773, 2264.482910156,  -20.000000000 )"
  	echo zoom $obj, 50
	echo png $fileDir/$obj"_spec.png"

	echo hide everything
	echo "delete axis*"
	echo load $PWD/src/PyMolUtils/LineUtils.py
	echo show cgo

	echo align "'"$relObj"'", "'"$obj"'"
	echo "hide (hydro)"
	echo "set solvent_radius, 4   
			alter all, vdw=6 
			sort
			set surface_quality, 1"
	echo show spheres
	echo color forest, $obj
	echo color teal, $relObj
	# echo set cartoon_cylindrical_helices, 1
	# echo set cartoon_transparency, 0.4, $relObj
	echo set sphere_transparency, 0.7

	# for chars
	# echo "set_view (\
	#      0.792623580,    0.107972473,   -0.600072622,\
	#      0.000102258,    0.984166384,    0.177220225,\
	#      0.609710038,   -0.140525863,    0.780067921,\
	#      0.000217728,   -0.001177430, -1687.383666992,\
	#     -6.746949673,   93.794700623,  -40.651931763,\
	#   	1256.293457031, 2118.493164062,  -20.000000000 )"
	# For weird shapes
	echo "set_view (\
		 0.682993114,   -0.467070848,    0.561571896,\
		 0.001815517,    0.769914448,    0.638143957,\
		-0.730422139,   -0.434827507,    0.526694298,\
		 0.000000000,    0.000000000, -1833.382080078,\
		 4.250000000,    6.750000000,    0.000000000,\
		1402.282592773, 2264.482910156,  -20.000000000 )"
  	echo zoom $obj, 50
	echo enable $obj
	echo png $fileDir/$obj"_design.png"
	echo enable $relObj
	echo disable $obj
	echo png $fileDir/$obj"_relaxed.png"
	echo
done