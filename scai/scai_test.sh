###
 # @file scai_test.sh
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
 # @brief ToDo: Missing description in ./scai_test.sh
 # @author Lauretta Schubert
 # @date 22.03.2016
###

# Set some envirionment variables to disable warning/error message

export SCAI_LOG=ERROR
export SCAI_TRACE=NONE
export SCAI_UNSUPPORTED=IGNORE

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
./common/test/commonTest

# Logging tests
(
    cd logging/test/
    echo "### loggingTest"
    ./test.sh
)

# Tracing tests
(
    cd tracing/test/
    echo "### tracingTest"
    ./test.sh
)

for CTX in ${CONTEXTS[*]}
do
	# Tasking tests
	echo "### taskingTest on ${CTX}"
	./tasking/test/taskingTest --SCAI_CONTEXT=${CTX}

	# HMemo tests
	echo "### hmemoTest on ${CTX}"
	./hmemo/test/hmemoTest --SCAI_CONTEXT=${CTX}

	# KRegistry tests
	echo "### kregistryTest on ${CTX}"
	./kregistry/test/kregistryTest --SCAI_CONTEXT=${CTX}

	# BLASKernel tests
	echo "### blaskernelTest on ${CTX}"
	./blaskernel/test/blaskernelTest --SCAI_CONTEXT=${CTX}

	# UtilsKernel tests
	echo "### utilskernelTest on ${CTX}"
	./utilskernel/test/utilskernelTest --SCAI_CONTEXT=${CTX}

	# SparseKernel tests
	echo "### sparsekernelTest on ${CTX}"
	./sparsekernel/test/sparsekernelTest --SCAI_CONTEXT=${CTX}

	# DMemo tests
	export SCAI_COMMUNICATOR=NO
	echo "### dmemoTest on ${CTX}"
	./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX}
	if [ "${MPI_FOUND}" != "" ]
	then
		export SCAI_COMMUNICATOR=MPI
		mpirun -np 1 ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX}
		mpirun -np 2 ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX}
		mpirun -np 3 ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX}
		mpirun -np 4 ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX}
	fi

	# LAMA tests
	echo "### lama_test on ${CTX}"
	( # lamaTest
		./lama/test/lamaTest --SCAI_CONTEXT=${CTX}
	)
	( # Storage Test
		./lama/test/storage/lamaStorageTest --SCAI_CONTEXT=${CTX}
	)
# # LAMA Distributed test removed?
#	(
#	    export SCAI_COMMUNICATOR=NO
#	    ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDistTest_${CTX}.xml --SCAI_CONTEXT=${CTX}
#	    if [ "${MPI_FOUND}" != "" ]
#	    then
#	        export SCAI_COMMUNICATOR=MPI
#	        mpirun -np 1 ./lama/test/distributed/lamaDistTest --SCAI_CONTEXT=${CTX}
#	        mpirun -np 2 ./lama/test/distributed/lamaDistTest --SCAI_CONTEXT=${CTX}
#	        mpirun -np 3 ./lama/test/distributed/lamaDistTest --SCAI_CONTEXT=${CTX}
#	        mpirun -np 4 ./lama/test/distributed/lamaDistTest --SCAI_CONTEXT=${CTX}
#	    fi
#    )
 	( # Matrix Test   
	    export SCAI_COMMUNICATOR=NO
	    ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX}
	    if [ "${MPI_FOUND}" != "" ]
	    then
	        export SCAI_COMMUNICATOR=MPI
	        mpirun -np 1 ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX}
	        mpirun -np 2 ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX}
	        mpirun -np 3 ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX}
	        mpirun -np 4 ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX}
	    fi
    )

	# Solver tests

	echo "### solverTest on ${CTX}"
	export SCAI_COMMUNICATOR=NO
	./solver/test/solverTest --SCAI_CONTEXT=${CTX}
	export SCAI_COMMUNICATOR=NO
	./solver/test/distributed/solverDistTest --SCAI_CONTEXT=${CTX}
	if [ "${MPI_FOUND}" != "" ]
	then
		export SCAI_COMMUNICATOR=MPI
		mpirun -np 1 ./solver/test/distributed/solverDistTest --SCAI_CONTEXT=${CTX}
		mpirun -np 2 ./solver/test/distributed/solverDistTest --SCAI_CONTEXT=${CTX}
		mpirun -np 3 ./solver/test/distributed/solverDistTest --SCAI_CONTEXT=${CTX}
		mpirun -np 4 ./solver/test/distributed/solverDistTest --SCAI_CONTEXT=${CTX}
	fi
done

unset CONTEXTS
