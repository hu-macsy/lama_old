###
 # @file lama/examples/run_all.sh
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief ToDo: Missing description in ./lama/examples/run_all.sh
 # @author Jan Ecker
 # @date 10.02.2016
###

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
RUN 1 bench/matvecmul.exe
RUN 1 bench/maxnorm.exe
RUN 1 bench/rowcol.exe
RUN 1 bench/sort.exe 10000
RUN 1 bench/scan.exe 10000

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
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi


#reset i
i=0

# run examples tutorial/*
RUN 1 tutorial/blas1.exe
RUN 1 tutorial/matrix.exe
RUN 1 tutorial/matrix1.exe
RUN 1 tutorial/scalar.exe
RUN 1 tutorial/simple.exe
RUN 1 tutorial/vector.exe
RUN 1  tutorial/vector_exp.exe

# check if there are unkown examples
count=`ls -l -la $MYDIR/tutorial/*.exe | wc -l`
if [ $count -ne $i ]; then
    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
    exit 1
fi
