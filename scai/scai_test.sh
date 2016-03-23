#!/bin/bash -e
# This script runs all test of the SCAI projects

# Set some envirionment variables to disable warning/error message

export SCAI_LOG=ERROR
export SCAI_TRACE=NONE
export SCAI_UNSUPPORTED=IGNORE

MPI_FOUND=$(which mpirun > /dev/null 2> /dev/null)

# Common tests

echo "### commonTest"
./common/test/commonTest

if [ -d cuda ];
then
    echo "### commonCUDATest"
    ./common/test/cuda/commonCUDATest
fi

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

# Tasking tests

echo "### taskingTest"
./tasking/test/taskingTest

# HMemo tests

echo "### hmemoTest"
./hmemo/test/hmemoTest

if [ -d cuda ];
then
    echo "### hmemoCUDATest"
    ./hmemo/test/cuda/hmemoCUDATest
fi

if [ -d mic ];
then
    echo "### hmemoMICTest"
    ./hmemo/test/mic/hmemoMICTest
fi

# KRegistry tests

echo "### kregistryTest"
./kregistry/test/kregistryTest

# BLASKernel tests

echo "### blaskernelTest"
./blaskernel/test/blaskernelTest

# UtilsKernel tests

echo "### utilskernelTest"
./utilskernel/test/utilskernelTest

# SparseKernel tests

echo "### sparsekernelTest"
./sparsekernel/test/sparsekernelTest

# DMemo tests

export SCAI_COMMUNICATOR=NO
echo "### dmemoTest"
./dmemo/test/dmemoTest
if [ "${MPI_FOUND}" != "" ]
then
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./dmemo/test/dmemoTest
    mpirun -np 2 ./dmemo/test/dmemoTest
    mpirun -np 3 ./dmemo/test/dmemoTest
    mpirun -np 4 ./dmemo/test/dmemoTest
fi

# LAMA tests

echo "### lama_test"
./lama/test/lamaTest
# TODO: should be removed
if [ -d distributed ]
then
    export SCAI_COMMUNICATOR=NO
    ./lama/test/distributed/lamaDistTest
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./lama/test/distributed/lamaDistTest
        mpirun -np 2 ./lama/test/distributed/lamaDistTest
        mpirun -np 3 ./lama/test/distributed/lamaDistTest
        mpirun -np 4 ./lama/test/distributed/lamaDistTest
    fi
fi

# Solver tests

echo "### solverTest"
export SCAI_COMMUNICATOR=NO
./solver/test/solverTest
export SCAI_COMMUNICATOR=NO
./solver/test/distributed/solverDistTest
if [ "${MPI_FOUND}" != "" ]
then
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./solver/test/distributed/solverDistTest
    mpirun -np 2 ./solver/test/distributed/solverDistTest
    mpirun -np 3 ./solver/test/distributed/solverDistTest
    mpirun -np 4 ./solver/test/distributed/solverDistTest
fi
