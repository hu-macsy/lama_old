#!/bin/bash
set -e

# clean up
make clean

# build examples
make

# run examples
./Aliasing.exe
./Array.exe
./BenchArray.exe
./BenchContext.exe
./Create.exe
./Threading.exe

# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne 6 ]; then
    echo "There are unkown executables (examples) in scai-hmemo, please add them to jenkins!"
    exit 1
fi

# run cuda tests
cd cuda

# clean up
make clean

# build examples
make

# run examples
#./AliasProblem.exe
./Allocate.exe
./CUBlasExample.exe
./CUDABenchContext.exe
./CUSparseExample.exe
./Devices.exe
./Example1.exe
#./Example2.exe
./MemBandwidth.exe
./Prefetch.exe

# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne 10 ]; then
    echo "There are unkown executables (examples) in scai-hmemo, please add them to jenkins!"
    exit 1
fi

cd ..