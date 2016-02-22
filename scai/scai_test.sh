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

# Common tests

(
    cd build/common/test
    ./CommonTest
)

# Logging tests

(
    cd build/logging/test
    ./test.sh
)

# HMemo tests

(
    cd build/hmemo/test
    ./MemoryTest
)

# KRegistry tests

(
	cd build/kregistry/test
	./KernelInterfaceTest
)

# BLASKernel tests

(
	cd build/blaskernel/test
	./BLASKernelTest
)

# UtilsKernel tests

(
	cd build/utilskernel/test
	./UtilsKernelTest
)

# SparseKernel tests

(
	cd build/sparsekernel/test
	./SparseKernelTest
)

# DMemo tests

(
    cd build/dmemo/test
    export SCAI_COMMUNICATOR=NO
    ./dmemoTest
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./dmemoTest
    mpirun -np 2 ./dmemoTest
    mpirun -np 3 ./dmemoTest
    mpirun -np 4 ./dmemoTest
)

# LAMA tests

(
    cd build/lama/test
    ./lama_test --log_level=test_suite | tee out.test
    cd distributed
    export SCAI_COMMUNICATOR=NO
    ./lama_dist_test
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./lama_dist_test
    mpirun -np 2 ./lama_dist_test
    mpirun -np 3 ./lama_dist_test
    mpirun -np 4 ./lama_dist_test
)

# Solver tests

( 
    cd build/solver/test
    export SCAI_COMMUNICATOR=NO
    ./SolverTest
    cd distributed
    export SCAI_COMMUNICATOR=NO
    ./SolverDistTest
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./SolverDistTest
    mpirun -np 2 ./SolverDistTest
    mpirun -np 3 ./SolverDistTest
    mpirun -np 4 ./SolverDistTest
)
