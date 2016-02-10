#!/bin/bash
set -e

# clean up
make clean
rm -rf *.exe *.o *.frm *.mtx *.frv *.amg *.vec *.trace

# build examples
make

# run examples

# bench/*
bench/conversion.exe
bench/matadd.exe
bench/matmul.exe

# check if there are unkown examples
count=`ls -l -la bench/*.exe | wc -l`
if [ $count -ne 3 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi

# labelrank/*
labelrank/labelrank.exe labelrank/affinity.mtx labelrank/labels.mtx

# check if there are unkown examples
count=`ls -l -la labelrank/*.exe | wc -l`
if [ $count -ne 1 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi

# tutorial/*
tutorial/blas1.exe
tutorial/matrix.exe
tutorial/scalar.exe
tutorial/simple.exe
tutorial/vector.exe
# currently disabled!!!!! (change check when enabling again!)
#tutorial/vector_exp.exe

# check if there are unkown examples
count=`ls -l -la tutorial/*.exe | wc -l`
if [ $count -ne 5 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi