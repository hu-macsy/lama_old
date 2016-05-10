###
 # @file scai_test_xml.sh
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
 # @brief ToDo: Missing description in ./scai_test_xml.sh
 # @author Lauretta Schubert
 # @date 22.03.2016
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

MPI_FOUND=$(which mpirun > /dev/null 2> /dev/null)

# Common tests

echo "### commonTest"
./common/test/commonTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/commonTest.xml

if [ -d cuda ];
then
    echo "### commonCUDATest"
    ./common/test/cuda/commonCUDATest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/commonCUDATest.xml
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
./tasking/test/taskingTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/taskingTest.xml

# HMemo tests

echo "### hmemoTest"
./hmemo/test/hmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoTest.xml

if [ -d cuda ];
then
    echo "### hmemoCUDATest"
    ./hmemo/test/cuda/hmemoCUDATest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoCUDATest.xml
fi

if [ -d mic ];
then
    echo "### hmemoMICTest"
    ./hmemo/test/mic/hmemoMICTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/hmemoMICTest.xml
fi

# KRegistry tests

echo "### kregistryTest"
./kregistry/test/kregistryTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/kregistryTest.xml

# BLASKernel tests

echo "### blaskernelTest"
./blaskernel/test/blaskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/kregistryTest.xml

# UtilsKernel tests

echo "### utilskernelTest"
./utilskernel/test/utilskernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/utilskernelTest.xml

# SparseKernel tests

echo "### sparsekernelTest"
./sparsekernel/test/sparsekernelTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/sparsekernelTest.xml

# DMemo tests

export SCAI_COMMUNICATOR=NO
echo "### dmemoTest"
./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemoTest.xml
if [ "${MPI_FOUND}" != "" ]
then
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo1Test.xml
    mpirun -np 2 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo2Test.xml
    mpirun -np 3 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo3Test.xml
    mpirun -np 4 ./dmemo/test/dmemoTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/dmemo4Test.xml
fi

# LAMA tests

echo "### lama_test"
./lama/test/lamaTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaTest.xml
# TODO: should be removed
if [ -d distributed ]
then
    export SCAI_COMMUNICATOR=NO
    ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDistTest.xml
    if [ "${MPI_FOUND}" != "" ]
    then
        export SCAI_COMMUNICATOR=MPI
        mpirun -np 1 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist1Test.xml
        mpirun -np 2 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist2Test.xml
        mpirun -np 3 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist3Test.xml
        mpirun -np 4 ./lama/test/distributed/lamaDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/lamaDist4Test.xml
    fi
fi

# Solver tests

echo "### solverTest"
export SCAI_COMMUNICATOR=NO
./solver/test/solverTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverTest.xml
export SCAI_COMMUNICATOR=NO
./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDistTest.xml
if [ "${MPI_FOUND}" != "" ]
then
    export SCAI_COMMUNICATOR=MPI
    mpirun -np 1 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist1Test.xml
    mpirun -np 2 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist2Test.xml
    mpirun -np 3 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist3Test.xml
    mpirun -np 4 ./solver/test/distributed/solverDistTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/solverDist4Test.xml
fi
