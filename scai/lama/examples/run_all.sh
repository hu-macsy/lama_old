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
echo "====================================================="
echo "==  Building and executing all scai lama examples  =="
echo "====================================================="
echo ""

cd $MYDIR

# build examples
make

# Use a counter to keep track of the number of executed examples
i=0

# run examples bench/*
RUN 1 bench/conversion.exe
RUN 1 bench/matadd.exe
RUN 1 bench/matmul.exe
RUN 1 bench/matnorm.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/bench/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi


#reset i
i=0

# run examples labelrank/*
RUN 1 labelrank/labelrank.exe $MYDIR/labelrank/affinity.mtx $MYDIR//labelrank/labels.mtx

# check if there are unkown examples
count=`ls -l -la $MYDIR/labelrank/*.exe | wc -l`
if [ $count -ne 1 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi


#reset i
i=0

# run examples tutorial/*
RUN 1 tutorial/blas1.exe
RUN 1 tutorial/matrix.exe
RUN 1 tutorial/scalar.exe
RUN 1 tutorial/simple.exe
RUN 1 tutorial/vector.exe
# currently disabled!!!!! (change check when enabling again!)
#RUN 1  tutorial/vector_exp.exe
i=$((i+1))


# check if there are unkown examples
count=`ls -l -la $MYDIR/tutorial/*.exe | wc -l`
if [ $count -ne 5 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi