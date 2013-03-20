#
#  @file cca.sh
# 
#  @license
#  Copyright (c) 2011
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
#  $
#

# Setting some enviroment variables

#export LAMA_LOG=config
export LAMA_UNSUPPORTED=IGNORE

echo `pwd`

path_lcov=/home/lama/bin/lcov
path_genhtml=/home/lama/bin/genhtml

# Creating dir named by YEARS_MONTHS_DAYS-HOURSMINUTES
dirname=$(date +%y_%m_%d-%k%M)
echo "Create coverage directory: ${dirname}"
mkdir ${dirname}

cd ${dirname}

# Clearing up environment
${path_lcov} --base-directory ../.. --directory ../.. --zerocounters



cd ..

# Running tests serial
echo "Running serial tests"
./lama_test

if [ -d distributed ];
then
	# Running parallel tests serial and with two processes
	cd distributed
	echo "Running distributed tests serial"
	./lama_dist_test
	echo "Running distributed tests with 2 processes"
	mpirun -np 2 ./lama_dist_test
	cd ..
fi

if [ -d cuda ];
then
	#Running CUDA tests
	cd cuda
	echo "Running cuda tests"
	./lama_cuda_test
	cd ..
fi

cd ${dirname}

#Running lcov and creating data
${path_lcov} --base-directory ../.. --directory ../.. --capture --output-file=data.info

#Extracting just Sourcefiles from LAMA/src/*
${path_lcov} --extract data.info "*/LAMA/src/*" --output-file=data.info

# Generating html-structure
${path_genhtml} data.info
