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

cd $MYDIR

# build examples
make

# Use a counter to keep track of the number of executed examples
i=0

# run examples
RUN 1 vectorAdd.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi

if [ -d $MYDIR/cuda ];
then

	cd $MYDIR/cuda

	# build examples
	make


	# reset counter for CUDA examples
	i=0

	# run CUDA examples
	#RUN 1 cuda/AliasProblem.exe

	# check if there are unkown examples
	count=`ls -l -la $MYDIR/cuda/*.exe | wc -l`
	if [ $count -ne $i ]; then
	    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
	    exit 1
	fi

	cd ..

fi

if [ -d $MYDIR/mic ];
then

	cd $MYDIR/mic

	# build examples
	make

	# reset counter for CUDA examples
	i=0

	# run mic examples
	#RUN 1 mic/MICExample1.exe

	# check if there are unkown examples
	count=`ls -l -la $MYDIR/mic/*.exe | wc -l`
	if [ $count -ne $i ]; then
	    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
	    exit 1
	fi

	cd ..
fi