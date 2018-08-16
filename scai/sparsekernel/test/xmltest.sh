###
 # @file xmltest.sh
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
 # @brief This file is a shellscript, which executes all sparsekernel tests and creates
 #        xml result files for further usage
 # @author Thomas Branes
 # @date 08.05.2013
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE

dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite


# Test for Host context is mandatory, other contexts are optional and tested if supported

SCAI_CONTEXT_LIST="Host"

if [ -d ../cuda ];
then
    SCAI_CONTEXT_LIST="${SCAI_CONTEXT_LIST} CUDA"
fi


# Options specific for Boost Unit Test set for all test runs

BOOST_TEST_ARGS="--output_format=XML --log_level=test_suite --report_level=no"

for ctx in ${SCAI_CONTEXT_LIST};
do
    echo "Running sparsekernel tests on ${ctx} context (no MKL/cuSparse)"
    SCAI_USE_MKL=0 SCAI_CUDA_USE_CUSPARSE=0 ./sparsekernelTest --SCAI_CONTEXT=${ctx} ${BOOST_TEST_ARGS} 1> ${dirname}/SparseKernel${ctx}Test.xml

    echo "Running sparsekernel tests on ${ctx} context (with MKL/cuSparse)"
    SCAI_USE_MKL=1 SCAI_CUDA_USE_CUSPARSE=1 ./sparsekernelTest --SCAI_CONTEXT=${ctx} ${BOOST_TEST_ARGS} 1> ${dirname}/SparseKernel${ctx}Test.xml
done
