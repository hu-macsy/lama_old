###
 # @file xmltest.sh
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
 # @brief This file is a shellscript, which executes all lama tests and creates
 #        xml result files for further usage
 # @author Thomas Brandes, Jan Ecker
 # @date 08.05.2013
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE

dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

# For configuration tests do not stop this script in case of error

set +e

# Test for Host context is mandatory, other contexts are optional and tested if supported

SCAI_CONTEXT_LIST="Host"
SCAI_OTHER_CONTEXTS="CUDA"

# check which of the optional contexts are supported in the current configuration

for ctx in ${SCAI_OTHER_CONTEXTS};
do
    # just jun one very simple test to see if context is supported

    ./lamaTest --SCAI_CONTEXT=${ctx} --run_test=VersionTest > /dev/null >& /dev/null

    if [ "$?" -eq 0 ]
    then
       SCAI_CONTEXT_LIST="${SCAI_CONTEXT_LIST} ${ctx}" 
    fi
done

# check if installed LAMA version supports MPI communicator

./lamaTest --SCAI_COMMUNICATOR=MPI --run_test=VersionTest > /dev/null >& /dev/null

if [ "$?" -eq 0 ]
then
   MPI_PROCS="2 3 4"
else
   MPI_PROCS=""
fi

# Options specific for Boost Unit Test set for all test runs

BOOST_TEST_ARGS="--output_format=XML --log_level=test_suite --report_level=no"

# For running unit tests do stop this script in case of an error

set -e 

echo "LAMA unit tests: MPI_PROCS=${MPI_PROCS}, SCAI_CONTEXT_LIST=${SCAI_CONTEXT_LIST}"

#Running tests for each supported context

for ctx in ${SCAI_CONTEXT_LIST};
do
    echo "Running lama tests on ${ctx} context"
    ./lamaTest --SCAI_CONTEXT=${ctx} ${BOOST_TEST_ARGS} 1> ${dirname}/lama${ctx}Test.xml

    echo "Running lama storage tests on ${ctx} context"
    ./storage/lamaStorageTest --SCAI_CONTEXT=${ctx} ${BOOST_TEST_ARGS} 1> ${dirname}/lamaStorage${ctx}Test.xml

    echo "Running lama matrix tests on ${ctx} context"
    ./matrix/lamaMatrixTest --SCAI_CONTEXT=${ctx} ${BOOST_TEST_ARGS} 1> ${dirname}/lamaMatrix${ctx}Test.xml

done

# Running parallel tests with different number of processes, but with only one thread per process

for i in ${MPI_PROCS} ;
do 
    echo "Running lama tests with $i processes"
    mpirun -np $i --output-filename ${dirname}/lamaMPITest.xml ./lamaTest --SCAI_COMMUNICATOR=MPI --SCAI_NUM_THREADS=1 ${BOOST_TEST_ARGS}
    echo "Running lama matrix tests with $i processes"
    mpirun -np $i --output-filename ${dirname}/lamaMatrixMPITest.xml ./matrix/lamaMatrixTest --SCAI_COMMUNICATOR=MPI --SCAI_NUM_THREADS=1 ${BOOST_TEST_ARGS}
done

echo "Tests finished, results in directory: ${dirname}"
