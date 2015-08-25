#!/bin/bash

# Skript to build all SCAI projects

SCAI_ROOT=/home/lschubert/localInstall/lamaStandalone

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
        # echo "cmake .. -DCMAKE_INSTALL_PREFIX=${SCAI_ROOT}\""
        cmake .. -DCMAKE_INSTALL_PREFIX=${SCAI_ROOT}
        make -j 4

        make install
        cd ../..
    done
fi

