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

# Options specific for Boost Unit Test set for all test runs

BOOST_TEST_ARGS="--output_format=XML --log_level=test_suite --report_level=no"

# Running dmemo tests (only Host)
echo "Running dmemo tests with NoCommunicator"
./dmemoTest --SCAI_COMMUNICATOR=NO ${BOOST_TEST_ARGS} 1> ${dirname}/dmemoTest.xml

if [ -d ../mpi ];
then
    echo "Running dmemo tests distributed with MPI"
	for i in 1 2 3 4 6 11;
	do
    	echo "Running MPI tests with $i processes"
    	mpirun -np $i --output-filename ${dirname}/dist_tests_mpi_${i}.xml ./dmemoTest --SCAI_COMMUNICATOR=MPI ${BOOST_TEST_ARGS} --SCAI_NUM_THREADS=1
    done
fi
