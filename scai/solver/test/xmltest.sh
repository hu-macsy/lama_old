###
 # @file xmltest.sh
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
./solverTest --SCAI_CONTEXT=CUDA --run_test=VersionTest > /dev/null >& /dev/null
if [ "$?" -eq 0 ]
then
   LAMA_SUPPORTS_CUDA=1
fi

# check if installed LAMA version supports MPI communicator

LAMA_SUPPORTS_MPI=0
./solverTest --SCAI_COMMUNICATOR=MPI --run_test=VersionTest > /dev/null >& /dev/null
if [ "$?" -eq 0 ]
then
   LAMA_SUPPORTS_MPI=1
fi

echo "LAMA SUPPORTS MPI = ${LAMA_SUPPORTS_MPI}, CUDA = ${LAMA_SUPPORTS_CUDA}"

# Stop execution as soon as test fails

set +e

BOOST_TEST_ARGS="--output_format=XML --log_level=${ERROR_LEVEL} --report_level=no"

# Running solver test (only Host)
echo "Running solver tests on Host"
./solverTest ${BOOST_TEST_ARGS} 1> ${dirname}/solverTest.xml

if [ "${LAMA_SUPPORTS_MPI}" -eq 1 ];
then
	for i in 2 3 4;
	do
        LAMA_ARGS="--SCAI_NUM_THREADS=1 --SCAI_CONTEXT=Host --SCAI_COMMUNICATOR=MPI"
    	echo "Running solver tests with $i processes"
        mpirun -np $i --output-filename ${dirname}/solverTest_mpi_$i.xml ./solverTest ${LAMA_ARGS} ${BOOST_TEST_ARGS}
    done
fi

echo "All tests are finished, results in directory: ${dirname}"
