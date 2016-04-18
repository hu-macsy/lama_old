#
#  @file code_coverage.sh
# 
#  @license
#  Copyright (c) 2009-2013
#  Fraunhofer Institute for Algorithms and Scientific Computing SCAI
#  for Fraunhofer-Gesellschaft
# 
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#  @endlicense
# 
#  @brief This file is a shellscript, which contains all necessary steps to 
#         measure code coverage of LAMA.
#  @author: Alexander Büchel, Lauretta Schubert
#  @date 15.08.2012
#  @since 1.0.0
#

#!/bin/bash -e

function create_dir
{
	# Creating dir using current unix timestamp
	local dirname=coverage_$(date +%s)
	mkdir ${dirname}
	echo "${dirname}"
}

# search for newest coverage_dir
function find_dir
{
	local dirs=($(ls -d coverage_*))

	biggest=0
	for _dir in ${dirs[@]}; do
		my_val=$(echo "${_dir}" | cut -d'_' -f2)
		if [ ${my_val} -ge ${biggest} ]; then
			biggest=${my_val}
		fi
	done
	echo "coverage_${biggest}"
}

# $1: code coveage dirname $2: base coverage dir
function prepare_coverage
{
	# Clearing up environment
	cd $1
	lcov --base-directory $2 --directory $2 --zerocounters
	cd -
}

# $1: code coveage dirname $2: base coverage dir $3: current source dir
function do_coverage
{
	cd $1

	local error_count=0

	#Running lcov and creating data
	lcov --base-directory $2 --directory $2 --capture --output-file=data.info
	if [ $? -ne 0 ]
	then
		echo "ERROR in running lcov"
		error_count=$(($error_count + $?))
	fi

	#Extracting just Sourcefiles
	lcov --extract data.info "$3/*" --output-file=data.info
	if [ $? -ne 0 ]
	then
		echo "ERROR in extracting lcov data"
		error_count=$(($error_count + $?))
	fi

	# Generating html-structure
	genhtml data.info
	if [ $? -ne 0 ]
	then
		echo "ERROR in generating html-structure"
		error_count=$(($error_count + $?))
	fi

	if [[ $error_count != 0 ]]
	then
	    exit 1
	fi

	cd -
}
