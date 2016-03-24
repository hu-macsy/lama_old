#
#  @file xmltest.sh
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
#  @brief This file is a shellscript, which executes all lama tests and creates
#         xml result files for further usage
#  @author: Jan Ecker
#  @date 08.05.2013
#

#!/bin/bash

# Creating dir named by YEAR_MONTH_DAY-HOURMINUTE
dirname=xmlresult_$(date +%s)
echo "Create result directory: ${dirname}"
mkdir ${dirname}

ERROR_LEVEL=test_suite

# Running kregistry serial (only Host)
echo "Running kregistry tests on Host"
./kregistryTest --output_format=XML --log_level=${ERROR_LEVEL} --report_level=no 1>${dirname}/KernelInterfaceTest.xml
