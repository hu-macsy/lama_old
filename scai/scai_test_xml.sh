###
 # @file scai_test_xml.sh
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
 # @brief Script to run all tests with XML result files
 # @author Lauretta Schubert
 # @date 22.03.2016
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

declare -a CONTEXTS

if [ -d dmemo/mpi ]
then
    # Lama has been built with MPI support so it should be available
    MPI_PROCS="2 4"
else
    MPI_PROCS=""
fi

# Use Host Context
CONTEXTS+=(Host)

# Use CUDA
if [ -d hmemo/cuda ]
then
	CONTEXTS+=(CUDA)
fi

# Define some recurring arguments

BOOST_TEST_ARGS="--output_format=XML --log_level=${ERROR_LEVEL} --report_level=no"
MPI_TEST_ARGS="--SCAI_COMMUNICATOR=MPI --SCAI_NUM_THREADS=1"

echo "## Used Contexts: ${CONTEXTS[*]}, MPI_PROCS=${MPI_PROCS}"

# Common tests
echo "### commonTest"
./common/test/commonTest ${BOOST_TEST_ARGS} 1> ${dirname}/commonTest.xml

# Logging tests
echo "### loggingTest"
(
    cd logging/test/
    ./test.sh
)

# Tracing tests
echo "### tracingTest"
(
    cd tracing/test/
    ./test.sh
)
# KRegistry tests
echo "### kregistryTest"
(
    ./kregistry/test/kregistryTest ${BOOST_TEST_ARGS} 1> ${dirname}/kregistryTest.xml
)
	
for CTX in ${CONTEXTS[*]}
do
    export SCAI_CONTEXT=${CTX}

	# Tasking tests
	echo "### taskingTest on ${CTX}"
	./tasking/test/taskingTest ${BOOST_TEST_ARGS} 1> ${dirname}/taskingTest_${CTX}.xml
	
	# HMemo tests
	echo "### hmemoTest on ${CTX}"
	./hmemo/test/hmemoTest ${BOOST_TEST_ARGS} 1> ${dirname}/hmemoTest_${CTX}.xml 
	
	# BLASKernel tests
	echo "### blaskernelTest on ${CTX}"
	./blaskernel/test/blaskernelTest ${BOOST_TEST_ARGS} 1> ${dirname}/kregistryTest_${CTX}.xml
	
	# UtilsKernel tests
	echo "### utilskernelTest on ${CTX}"
	./utilskernel/test/utilskernelTest ${BOOST_TEST_ARGS} 1>${dirname}/utilskernelTest_${CTX}.xml
	
	# SparseKernel tests
	echo "### sparsekernelTest on ${CTX}"
	./sparsekernel/test/sparsekernelTest ${BOOST_TEST_ARGS} 1>${dirname}/sparsekernelTest_${CTX}.xml
	
	# DMemo tests
	echo "### dmemoTest on ${CTX}"
	./dmemo/test/dmemoTest ${BOOST_TEST_ARGS} --SCAI_COMMUNICATOR=NO 1>${dirname}/dmemoTest_${CTX}.xml
    for np in ${MPI_PROCS} ;
    do 
	    echo "### dmemoTest on ${CTX} with ${np} MPI processes"
	    mpirun -np ${np} -output-filename ${dirname}/dmemoTest_${CTX}_${np}.xml ./dmemo/test/dmemoTest ${BOOST_TEST_ARGS} ${MPI_TEST_ARGS} 
    done
	
	# LAMA tests
	echo "### lamaTest on ${CTX}"
	( # lamaTest
		./lama/test/lamaTest ${BOOST_TEST_ARGS} 1>${dirname}/lamaTest_${CTX}.xml
        for np in ${MPI_PROCS} ;
        do 
	        echo "### lamaTest on ${CTX} with ${np} MPI processes"
	        mpirun -np ${np} -output-filename ${dirname}/lamaTest_${CTX}_${np}.xml ./lama/test/lamaTest ${BOOST_TEST_ARGS} ${MPI_TEST_ARGS} 
        done
	)
	echo "### lamaStorageTest on ${CTX}"
	( # Storage Test
		./lama/test/storage/lamaStorageTest ${BOOST_TEST_ARGS} --SCAI_COMMUNICATOR=NO 1> ${dirname}/lamaStorageTest_${CTX}.xml
	)
	echo "### lamaMatrixTest on ${CTX}"
 	( # Matrix Test   
	    ./lama/test/matrix/lamaMatrixTest ${BOOST_TEST_ARGS} --SCAI_COMMUNICATOR=NO 1> ${dirname}/lamaMatrixTest_${CTX}.xml
        for np in ${MPI_PROCS} ;
        do 
	        echo "### lamaMatrixTest on ${CTX} with ${np} MPI processes"
	        mpirun -np ${np} -output-filename ${dirname}/lamaMatrixTest_${CTX}_${np}.xml ./lama/test/matrix/lamaMatrixTest ${BOOST_TEST_ARGS} ${MPI_TEST_ARGS} 
        done
    )
	
	# Solver tests
	echo "### solverTest on ${CTX}"
	./solver/test/solverTest ${BOOST_TEST_ARGS} --SCAI_COMMUNICATOR=NO 1> ${dirname}/solverTest_${CTX}.xml
    for np in ${MPI_PROCS} ;
    do 
        if [ "$CTX" == "CUDA" ]; then
            echo "skip MPI solver test with CUDA"
        else
	        echo "### solverTest on ${CTX} with ${np} MPI processes"
	        echo "mpirun -np ${np} ./solver/test/solverTest ${BOOST_TEST_ARGS} ${MPI_TEST_ARGS}"
	        mpirun -np ${np} -output-filename ${dirname}/solverTest_${CTX}_${np}.xml ./solver/test/solverTest ${BOOST_TEST_ARGS} ${MPI_TEST_ARGS} 
        fi
	done
done

echo "All tests are finished, for XML result files see ${dirname}"

if [ -n "${MPI_PROCS}" ]; then

	# rename the MPI output file to end with xml

	for filename in ${dirname}/*.xml.* ;
	do
	    # rename <filename>.xml.<rank> to <filename>.<rank>.xml
	    newname=${filename/\.xml/}.xml
	    mv $filename $newname
	done

fi

unset CONTEXTS
