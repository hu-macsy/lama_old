#!/bin/bash
set -e

# clean up
make clean

# build examples
make

# run examples
SCAI_TRACE=time ./TraceCalls.exe
SCAI_TRACE=time ./TraceTest.exe
SCAI_TRACE=time ./TraceThreads.exe

# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne 3 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi