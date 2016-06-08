#!/bin/bash -e

###
 # @file scai_code_coverage_functions.sh
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 # @brief This file is a shellscript, which contains all necessary steps to 
 #        measure code coverage of LAMA.
 # @author Alexander Büchel, Lauretta Schubert
 # @date 15.08.2012
###

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
