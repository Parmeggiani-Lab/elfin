To compile:
	make

To compile with your choice of compiler e.g. clang++:
	make CXX='clang++'

To get command line help:
	./bin/elfin -h

To run:
	./bin/elfin [arguments to override config.json]

To edit config:
	<your editor> config.json

Notes:
	1. Makefile adds jutil/src/ to include path, 
		that's why util.h is included as if it was in 
		the same directory everywhere

	2. If jutil/src is not properly installed, you need to:
		2.1) Go to root directory of elfin
		2.2) do "git submodule init; git submodule update"
