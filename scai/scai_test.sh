#!/bin/bash

# This script runs all test of the SCAI projects

# Skript to build all SCAI projects

function checkErrorValue( ) {
	$*
	if [ "$?" -ne 0 ];
	then
		echo "Build Failed. Aborting..."
		exit 1
	fi
}

# Set some envirionment variables to disable warning/error message

export SCAI_LOG=ERROR
export SCAI_TRACE=NONE
export SCAI_UNSUPPORTED=IGNORE

MPI_FOUND=$(which mpirun > /dev/null 2> /dev/null)

# Common tests

(
    cd common/test
    echo "### commonTest"
    ./commonTest

    if [ -d cuda ];
    then
        cd cuda
        echo "### commonCUDATest"
        ./commonCUDATest
    fi
)

# Logging tests

(
    cd logging/test
    echo "### loggingTest"
    ./test.sh
)

# Tracing tests

(
    cd tracing/test
    echo "### tracingTest"
    ./test.sh
)

# Tasking tests

(
    cd tasking/test
    echo "### taskingTest"
    ./taskingTest
)

# HMemo tests

(
    cd hmemo/test
    echo "### hmemoTest"
    ./hmemoTest

    (
    if [ -d cuda ];
    then
        cd cuda
        echo "### hmemoCUDATest"
        ./hmemoCUDATest
    fi
    )
    (
    if [ -d mic ];
    then
        cd mic
        echo "### hmemoMICTest"
        ./hmemoMICTest
    fi
    )
)

# KRegistry tests

(
	cd kregistry/test
    echo "### kregistryTest"
	./kregistryTest
)

# BLASKernel tests

(
	cd blaskernel/test
    echo "### blaskernelTest"
	./blaskernelTest
)

# UtilsKernel tests

(
	cd utilskernel/test
    echo "### utilskernelTest"
	./utilskernelTest
)

# SparseKernel tests

(
	cd sparsekernel/test
    echo "### sparsekernelTest"
	./sparsekernelTest
)

# DMemo tests

(
    cd dmemo/test
    export SCAI_COMMUNICATOR=NO
    echo "### dmemoTest"
    ./dmemoTest
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./dmemoTest
        mpirun -np 2 ./dmemoTest
        mpirun -np 3 ./dmemoTest
        mpirun -np 4 ./dmemoTest
    fi
)

# LAMA tests

(
    cd lama/test
    echo "### lama_test"
    ./lamaTest --log_level=test_suite | tee out.test
    // TODO: should be removed
    if [ -d distributed ]
    then
        cd distributed
        export SCAI_COMMUNICATOR=NO
        ./lamaDistTest
        if [ "${MPI_FOUND}" != "" ]
        then
            export SCAI_COMMUNICATOR=MPI
            mpirun -np 1 ./lamaDistTest
            mpirun -np 2 ./lamaDistTest
            mpirun -np 3 ./lamaDistTest
            mpirun -np 4 ./lamaDistTest
        fi
    fi
)

# Solver tests

( 
    cd solver/test
    echo "### solverTest"
    export SCAI_COMMUNICATOR=NO
    ./solverTest
    cd distributed
    export SCAI_COMMUNICATOR=NO
    ./solverDistTest
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./solverDistTest
        mpirun -np 2 ./solverDistTest
        mpirun -np 3 ./solverDistTest
        mpirun -np 4 ./solverDistTest
    fi
)
