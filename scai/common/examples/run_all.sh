#!/bin/bash
set -e

# clean up
make clean

# build examples
make

function checkErrorValue( ) {
    $*
    if [ "$?" -ne 0 ];
    then
        echo "Example: $cmd Failed. Aborting..."
        exit 1
    fi
}

# run examples

checkErrorValue ./Barrier.exe
checkErrorValue ./CriticalRegion.exe
checkErrorValue ./DemoComplex.exe
checkErrorValue ./DemoFactory.exe
checkErrorValue ./DemoFactory1.exe
checkErrorValue ./DemoFunction.exe
checkErrorValue ./DemoMath.exe
checkErrorValue ./DemoPointer.exe
checkErrorValue ./DemoSettings.exe
checkErrorValue ./DemoTypeTrait.exe
checkErrorValue ./ExceptionDemo.exe
checkErrorValue ./BenchPointers.exe 
checkErrorValue ./TimePrecision.exe 
checkErrorValue ./UseModule.exe ./DummyModule.so 
checkErrorValue ./UseModule.exe ./Module.so 

