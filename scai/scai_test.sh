###
 # @file scai_test.sh
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
 # @brief Script that runs all tests of all SCAI Projects
 # @author Lauretta Schubert
 # @date 22.03.2016
###

# Set some envirionment variables to disable warning/error message

export SCAI_LOG=ERROR
export SCAI_TRACE=NONE
export SCAI_UNSUPPORTED=IGNORE

if [ -d dmemo/mpi ]
then 
    MPI_PROCS="2 4"
else
    MPI_PROCS=""
fi

declare -a CONTEXTS

# Use Host Context
CONTEXTS+=(Host)

# Use CUDA
if [ -d hmemo/cuda ]
then
    CONTEXTS+=(CUDA)
fi

echo "## Used Contexts: ${CONTEXTS[*]}, MPI_PROCS=${MPI_PROCS}"

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

# KRegistry tests
(
    echo "### kregistryTest"
    ./kregistry/test/kregistryTest
)

for CTX in ${CONTEXTS[*]}
do
    # Tasking tests
    echo "### taskingTest on ${CTX}"
    ./tasking/test/taskingTest --SCAI_CONTEXT=${CTX}

    # HMemo tests
    echo "### hmemoTest on ${CTX}"
    ./hmemo/test/hmemoTest --SCAI_CONTEXT=${CTX}

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
    echo "### dmemoTest on ${CTX}"
    ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=NO
    for np in ${MPI_PROCS} ;
    do
        echo "### dmemoTest on ${CTX} with ${np} MPI processes"
        mpirun -np ${np} ./dmemo/test/dmemoTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=MPI
    done

    # LAMA tests
    echo "### lamaTest on ${CTX}"
    ( # lamaTest
        ./lama/test/lamaTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=NO
        for np in ${MPI_PROCS} ;
        do
            echo "### lamaTest on ${CTX} with ${np} MPI processes"
            mpirun -np ${np}  ./lama/test/lamaTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=MPI --SCAI_NUM_THREADS=1
        done
    )
    echo "### lamaStorageTest on ${CTX}"
    ( # Storage Test
        ./lama/test/storage/lamaStorageTest --SCAI_CONTEXT=${CTX}
    )
    echo "### lamaMatrixTest on ${CTX}"
    ( # Matrix Test   
        ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=NO
        for np in ${MPI_PROCS} ;
        do
            echo "### lamaMatrixTest on ${CTX} with ${np} MPI processes"
            mpirun -np ${np} ./lama/test/matrix/lamaMatrixTest --SCAI_CONTEXT=${CTX} --SCAI_NUM_THREADS=1 --SCAI_COMMUNICATOR=MPI
        done
    )

    # Solver tests

    echo "### solverTest on ${CTX}"
    ./solver/test/solverTest --SCAI_CONTEXT=${CTX} --SCAI_COMMUNICATOR=NO
    for np in ${MPI_PROCS} ;
    do
       echo "### solverTest on ${CTX} with ${np} MPI processes"
       mpirun -np ${np} ./solver/test/solverTest --SCAI_COMMUNICATOR=MPI --SCAI_NUM_THREADS=1
    done
done

unset CONTEXTS
