###
 # @file hmemo/examples/run_all.sh
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief ToDo: Missing description in ./hmemo/examples/run_all.sh
 # @author Jan Ecker
 # @date 10.02.2016
###

# Makes the bash exit if one commands returns with an error
set -e

# Get location of the script to properly call all example scripts
MYDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

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
RUN 1 Aliasing.exe
RUN 1 Array.exe
RUN 1 BenchArray.exe
RUN 1 BenchContext.exe
RUN 1 Create.exe
RUN 1 Threading.exe
RUN 1 RemoteAccess.exe
RUN 1 VectorArray.exe

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
	RUN 1 cuda/AliasProblem.exe
	RUN 1 cuda/Allocate.exe
	RUN 1 cuda/CUBlasExample.exe
	RUN 1 cuda/CUDABenchContext.exe
	RUN 1 cuda/CUSparseExample.exe
	RUN 1 cuda/Devices.exe
	RUN 1 cuda/Example1.exe
	RUN 1 cuda/Example2.exe
	RUN 1 cuda/MemBandwidth.exe
	RUN 1 cuda/Prefetch.exe

	# check if there are unkown examples
	count=`ls -l -la $MYDIR/cuda/*.exe | wc -l`
	if [ $count -ne $i ]; then
	    echo "There are unknown executables in this directory, please add all examples to the related run_all.sh script!"
	    exit 1
	fi

	cd ..

fi
