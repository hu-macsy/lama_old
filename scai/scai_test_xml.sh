###
 # @file scai_test_xml.sh
 #
 # @license
 # Copyright (c) 2009-2016
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
 # @endlicense
 #
 # @brief ToDo: Missing description in ./scai_test_xml.sh
 # @author Lauretta Schubert
 # @date 22.03.2016
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

# Use MPI
MPI_FOUND=$(which mpirun > /dev/null 2> /dev/null)

declare -a CONTEXTS

# Use Host Context
CONTEXTS+=(Host)

# Use CUDA
if [ -d hmemo/cuda ]
then
	CONTEXTS+=(CUDA)
fi

# Use MIC
if [ -d hmemo/mic ]
then
	CONTEXTS+=(MIC)
fi

echo "## Used Contexts: ${CONTEXTS[*]}"

# Common tests
echo "### commonTest"
./common/test/commonTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/commonTest.xml

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
	
for CTX in ${CONTEXTS[*]}
do
	# Tasking tests
	echo "### taskingTest on ${CTX}"
	./tasking/test/taskingTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/taskingTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# HMemo tests
	echo "### hmemoTest on ${CTX}"
	./hmemo/test/hmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# KRegistry tests
	echo "### kregistryTest on ${CTX}"
	./kregistry/test/kregistryTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/kregistryTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# BLASKernel tests
	echo "### blaskernelTest on ${CTX}"
	./blaskernel/test/blaskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/kregistryTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# UtilsKernel tests
	echo "### utilskernelTest on ${CTX}"
	./utilskernel/test/utilskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/utilskernelTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# SparseKernel tests
	echo "### sparsekernelTest on ${CTX}"
	./sparsekernel/test/sparsekernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/sparsekernelTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	
	# DMemo tests
	echo "### dmemoTest on ${CTX}"
	export SCAI_COMMUNICATOR=NO
	./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemoTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	if [ "${MPI_FOUND}" != "" ]
	then
	    export SCAI_COMMUNICATOR=MPI
	    mpirun -np 1 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo1Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 2 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo2Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 3 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo3Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 4 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo4Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	fi
	
	# LAMA tests
	echo "### lama_test on ${CTX}"
	( # lamaTest
		./lama/test/lamaTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	)
	( # Storage Test
		./lama/test/storage/lamaStorageTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaStorageTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	)
# # LAMA Distributed test removed?
#	(
#	    export SCAI_COMMUNICATOR=NO
#	    ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDistTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	    if [ "${MPI_FOUND}" != "" ]
#	    then
#	        export SCAI_COMMUNICATOR=MPI
#	        mpirun -np 1 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist1Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	        mpirun -np 2 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist2Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	        mpirun -np 3 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist3Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	        mpirun -np 4 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist4Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	    fi
#    )
 	( # Matrix Test   
	    export SCAI_COMMUNICATOR=NO
	    ./lama/test/matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrixTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    if [ "${MPI_FOUND}" != "" ]
	    then
	        export SCAI_COMMUNICATOR=MPI
	        mpirun -np 1 ./lama/test/matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrix1Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	        mpirun -np 2 ./lama/test/matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrix2Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	        mpirun -np 3 ./lama/test/matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrix3Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	        mpirun -np 4 ./lama/test/matrix/lamaMatrixTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaMatrix4Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    fi
    )
	
	# Solver tests
	echo "### solverTest on ${CTX}"
	export SCAI_COMMUNICATOR=NO
	./solver/test/solverTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	export SCAI_COMMUNICATOR=NO
	./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDistTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
	if [ "${MPI_FOUND}" != "" ]
	then
	    export SCAI_COMMUNICATOR=MPI
	    mpirun -np 1 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist1Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 2 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist2Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 3 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist3Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	    mpirun -np 4 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist4Test_${CTX}.xml --SCAI_CONTEXT=${CTX}
	fi
done


unset CONTEXTS
