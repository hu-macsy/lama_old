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
echo "======================================================="
echo "==  Building and executing all scai solver examples  =="
echo "======================================================="
echo ""

cd $MYDIR

# build examples
make

# cleanup
rm -rf lama.png
rm -rf gr_30_30.mtx.*.png
rm -rf *.exe *.o *.frm *.mtx *.frv *.amg *.vec *.trace

# Use a counter to keep track of the number of executed examples
i=0

# run examples lecture/*
RUN 1 lecture/task0.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task1a.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task1b.exe
RUN 1 lecture/task2.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task2a.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task3.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task4.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 lecture/task5.exe $MYDIR/lecture/gr_30_30.mtx

# check if there are unkown examples
count=`ls -l -la $MYDIR/lecture/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi



#reset counter
i=0

# run examples solver/*
cd $MYDIR/solver
RUN 1 solver/matrix_generator.exe example 3 27 100 100 100
RUN 1 solver/matrix_convert.exe -mm example.frv example.mtx
RUN 1 solver/vector_generator.exe example2.mtx 1000 1
RUN 1 solver/cg_solver.exe example
RUN 1 solver/gmres_solver.exe example
RUN 1 solver/amg_solver.exe example 3
RUN 1 solver/solver.exe example
RUN 0 solver/solver.exe example Jacobi 10
RUN 0 solver/solver.exe example GMRES 3
RUN 1 solver/lama_info.exe 

# check if there are unkown examples
count=`ls -l -la $MYDIR/solver/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi



#reset counter
i=0

# run examples spy/*
RUN 1 spy/amg_spy.exe $MYDIR/lecture/gr_30_30.mtx
RUN 1 spy/spy.exe $MYDIR/lecture/gr_30_30.mtx

# check if there are unkown examples
count=`ls -l -la $MYDIR/spy/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi


#reset counter
i=0

# run examples spy/*
RUN 1 myJacobi/myJacobi.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/myJacobi/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi
