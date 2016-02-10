#!/bin/bash
set -e

# clean up
make clean

# build examples
make

# run examples
./LogLevels.exe
./LogOpenMP.exe

# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne 2 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi