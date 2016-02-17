#!/bin/bash

# Makes the bash exit if one commands returns with an error
set -e

echo ""
echo "======================================================="
echo "==  Building and executing all scai common examples  =="
echo "======================================================="
echo ""

# clean up
make clean

# build examples
make

# Use a counter to keep track of the number of executed examples
i=0

# run examples
./Barrier.exe                       ;i=$((i+1))
./CriticalRegion.exe                ;i=$((i+1))
./DemoComplex.exe                   ;i=$((i+1))
./DemoFactory.exe                   ;i=$((i+1))
./DemoFactory1.exe                  ;i=$((i+1))
./DemoFunction.exe                  ;i=$((i+1))
./DemoMath.exe                      ;i=$((i+1))
./DemoPointer.exe                   ;i=$((i+1))
./DemoSettings.exe                  ;i=$((i+1))
./DemoTypeTrait.exe                 ;i=$((i+1))
./ExceptionDemo.exe                 ;i=$((i+1))
./BenchPointers.exe                 ;i=$((i+1))
./TimePrecision.exe                 ;i=$((i+1))
./UseModule.exe ./DummyModule.so    ;i=$((i+1))  # increment only once as they both use the UseModule.exe
./UseModule.exe ./Module.so         


# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi
