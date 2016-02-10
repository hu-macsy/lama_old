#!/bin/bash

# clean up
make clean

# build examples
make

# run examples
./Barrier.exe
./CriticalRegion.exe
./DemoComplex.exe
./DemoFactory.exe
./DemoFactory1.exe
./DemoFunction.exe
./DemoMath.exe
./DemoPointer.exe
./DemoSettings.exe
./DemoTypeTrait.exe
./ExceptionDemo.exe
./TimePrecision.exe
./UseModule.exe ./DummyModule.so
./UseModule.exe ./Module.so

# check if there are unkown examples
count=`ls -l -la *.exe | wc -l`
if [ $count -ne 13 ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi