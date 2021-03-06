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
 # @brief This file is a shellscript, which executes all tasking tests and creates
 #        xml result files for further usage
 # @author Thomas Brandes / Jan Ecker
 # @date 08.07.2016
###

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

# Running tasking tests (only Host)
echo "Running tasking tests on Host"
./taskingTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/taskingTest.xml

# Running tasking CUDA tests
if [ -d cuda ];
then
	echo "Running tasking tests for CUDA"
	./cuda/taskingCUDATest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/taskingCUDATest.xml
fi
