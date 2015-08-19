#!/bin/bash

# Skript to build all SCAI projects

SCAI_ROOT=/home/brandes/local/scai

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
        # echo "cmake .. -DCMAKE_INSTALL_PREFIX=${SCAI_ROOT} -DADDITIONAL_WARNING_FLAGS=\"-Wextra -Wall\""
        cmake .. -DCMAKE_INSTALL_PREFIX=${SCAI_ROOT}
        make install
        cd ../..
    done
fi

