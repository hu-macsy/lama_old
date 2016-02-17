#!/bin/bash

# Makes the bash exit if one commands returns with an error
set -e

# Get location of the script to properly call all example scripts
MYDIR="$(dirname "$(readlink -f "$0")")"

# Function that executes an example and count up a counter
# Usage: RUN COUNT[0|1] EXECUTABLE
#
function RUN ( ) {
    # count up for each new example
    i=$((i+$1))
    
    echo ""
    echo "Executing: ${@:2}"
    $MYDIR/${@:2}
}

echo ""
echo "======================================================"
echo "==  Building and executing all scai hmemo examples  =="
echo "======================================================"
echo ""

# build examples
make

# Use a counter to keep track of the number of executed examples
i=0

# run examples
RUN 1 Aliasing.exe
RUN 1 Array.exe
RUN 1 BenchArray.exe
RUN 1 BenchContext.exe
RUN 1 Create.exe
RUN 1 Threading.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi


# reset counter for CUDA examples
i=0

# run CUDA examples
#RUN 1 cuda/AliasProblem.exe
i=$((i+1))

RUN 1 cuda/Allocate.exe
RUN 1 cuda/CUBlasExample.exe
RUN 1 cuda/CUDABenchContext.exe
RUN 1 cuda/CUSparseExample.exe
RUN 1 cuda/Devices.exe
RUN 1 cuda/Example1.exe

#RUN 1 cuda/Example2.exe
i=$((i+1))

RUN 1 cuda/MemBandwidth.exe
RUN 1 cuda/Prefetch.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/cuda/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi