#!/bin/bash

# scai_build.sh clean
# scai_build.sh <install_dir>

# Skript to build all SCAI projects

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
        cmake .. -DCMAKE_INSTALL_PREFIX=$1 -DCMAKE_BUILD_TYPE=Release
        make -j 4

        make install
        cd ../..
    done
fi

