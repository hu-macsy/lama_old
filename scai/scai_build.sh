#!/bin/bash

# scai_build.sh clean
# scai_build.sh <install_dir>

# Skript to build all SCAI projects

function checkErrorValue( ) {
	$*
	if [ "$?" -ne 0 ];
	then
		echo "Build Failed. Aborting..."
		exit 1
	fi
}

if [ "$#" -ne 1 ]; then
   echo "Illegal call with $# arguments"
   echo "scai_build.sh clean"
   echo "scai_build.sh <install_prefix>"
   exit 1
fi

PROJECTS="common logging tracing tasking hmemo lama"

if [ "$1" == "clean" ]; then
    for project in $PROJECTS
    do
        echo "Clean SCAI project $project" 
        cd $project
        rm -rf build
        cd ..
    done
else
    for project in $PROJECTS
    do
        echo "Build SCAI project $project"
        cd $project
        mkdir -p build
        cd build
        # echo "cmake .. -DCMAKE_INSTALL_PREFIX=$1\""
        checkErrorValue cmake .. -DCMAKE_INSTALL_PREFIX=$1 -DCMAKE_BUILD_TYPE=Release 

	checkErrorValue make -j 4 
        checkErrorValue make install
        
	cd ../..
    done
fi

