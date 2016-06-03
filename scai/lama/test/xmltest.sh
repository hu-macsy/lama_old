#!/bin/bash -e

###
 # @file xmltest.sh
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief This file is a shellscript, which executes all lama tests and creates
 #        xml result files for further usage
 # @author Jan Ecker
 # @date 08.05.2013
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

set +e

# check if installed LAMA version supports CUDA context

LAMA_SUPPORTS_CUDA=0
./lamaTest --SCAI_CONTEXT=CUDA --run_test=VersionTest > /dev/null >& /dev/null
if [ "$?" -eq 0 ]
then
   LAMA_SUPPORTS_CUDA=1
fi

# check if installed LAMA version supports MPI communicator

LAMA_SUPPORTS_MPI=0
./lamaTest --SCAI_COMMUNICATOR=MPI --run_test=VersionTest > /dev/null >& /dev/null
if [ "$?" -eq 0 ]
then
   LAMA_SUPPORTS_MPI=1
fi

set -e 

echo "LAMA support checks: LAMA_SUPPORTS_MPI=${LAMA_SUPPORTS_MPI}, LAMA_SUPPORTS_CUDA=${LAMA_SUPPORTS_CUDA}"

# Running lama tests (only Host)
echo "Running lama tests on Host context"
./lamaTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaHostTest.xml

echo "Running lama storage tests on Host context"
./storage/lamaStorageTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaStorageHostTest.xml

echo "Running lama matrix tests on Host context"
./matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrixHostTest.xml

if [ "${LAMA_SUPPORTS_MPI}" -eq 1 ];
then
    # Running parallel tests serial and with two processes
    echo "Running matrix tests with 3 processes"
    mpirun -np 3 --output-filename ${dirname}/lamaMatrixMPITest.xml ./matrix/lamaMatrixTest --SCAI_COMMUNICATOR=MPI --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no

	#for i in 2 3 4;
	#do
    #	echo "Running distributed tests with $i processes"
    #	mpirun -np $i --output-filename ${dirname}/dist_tests_mpi.xml distributed/lamaDistTest --output_format=XML --log_level=all --report_level=no
    #done
fi

#Running CUDA tests

if [ "${LAMA_SUPPORTS_CUDA}" -eq 1 ];
then

    echo "Running lama tests on CUDA context"
    ./lamaTest --SCAI_CONTEXT=CUDA --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaCudaTest.xml

    echo "Running lama storage tests on CUDA context"
    ./storage/lamaStorageTest --SCAI_CONTEXT=CUDA --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaStorageCudaTest.xml

    echo "Running lama matrix tests on CUDA context"
    ./matrix/lamaMatrixTest --SCAI_CONTEXT=CUDA --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrixHostTest.xml
fi

