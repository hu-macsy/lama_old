#!/bin/bash -e
# This script runs all test of the SCAI projects

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite


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
        ./commonCUDATest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/commonTest.xml
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
    ./taskingTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/taskingTest.xml
)

# HMemo tests

(
    cd hmemo/test
    echo "### hmemoTest"
    ./hmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoTest.xml

    (
    if [ -d cuda ];
    then
        cd cuda
        echo "### hmemoCUDATest"
        ./hmemoCUDATest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoCUDATest.xml
    fi
    )
    (
    if [ -d mic ];
    then
        cd mic
        echo "### hmemoMICTest"
        ./hmemoMICTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoMICTest.xml
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
	./blaskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/kregistryTest.xml
)

# UtilsKernel tests

(
	cd utilskernel/test
    echo "### utilskernelTest"
	./utilskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/utilskernelTest.xml
)

# SparseKernel tests

(
	cd sparsekernel/test
    echo "### sparsekernelTest"
	./sparsekernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/sparsekernelTest.xml
)

# DMemo tests

(
    cd dmemo/test
    export SCAI_COMMUNICATOR=NO
    echo "### dmemoTest"
    ./dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemoTest.xml
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo1Test.xml
        mpirun -np 2 ./dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo2Test.xml
        mpirun -np 3 ./dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo3Test.xml
        mpirun -np 4 ./dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo4Test.xml
    fi
)

# LAMA tests

(
    cd lama/test
    echo "### lama_test"
    ./lamaTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaTest.xml
    # TODO: should be removed
    if [ -d distributed ]
    then
        cd distributed
        export SCAI_COMMUNICATOR=NO
        ./lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDistTest.xml
        if [ "${MPI_FOUND}" != "" ]
        then
            export SCAI_COMMUNICATOR=MPI
            mpirun -np 1 ./lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist1Test.xml
            mpirun -np 2 ./lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist2Test.xml
            mpirun -np 3 ./lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist3Test.xml
            mpirun -np 4 ./lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist4Test.xml
        fi
    fi
)

# Solver tests

( 
    cd solver/test
    echo "### solverTest"
    export SCAI_COMMUNICATOR=NO
    ./solverTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverTest.xml
    cd distributed
    export SCAI_COMMUNICATOR=NO
    ./solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDistTest.xml
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist1Test.xml
        mpirun -np 2 ./solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist2Test.xml
        mpirun -np 3 ./solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist3Test.xml
        mpirun -np 4 ./solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist4Test.xml
    fi
)
