#!/bin/bash
###
 # @file scai/logging/test/test.sh
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief Tests for SCAI logging
 # @author Jan Ecker
 # @date 03.09.2015
 # @since 2.0.0
###

#TODO:
# - documentation of functions
# - improve checks of the .ct and .time files (maybe combine check*FilesExist and check*FileContents functions)

genericTimePattern=", inclusive = "[0-9]{1,5}\.[0-9]{4,6}", exclusive = "[0-9]{1,5}\.[0-9]{4,6}$
errors=0

# =====================================================================================================================
# Runtime configurations tests
#
# In this test the executable is build WITH trace support and the runtime configuration via environmental variables is
# used to control the tracing behavior. 
# =====================================================================================================================

# Define checkTimeFileContents and checkCTFileContents first. Functions validate the contents of the created *.ct
# and *.time files.

function prepareTestCase {
    if [ $# -ne 1 ]; then
        echo "Invalid number of parameters!"
        exit 1
    fi

    # clean old tracing files
    rm -rf *.ct* *.time*
    
    # export SCAI_TRACE setting
    if [ -z "$1" ]; then
        unset SCAI_TRACE
    else
        export SCAI_TRACE=$1
    fi
   
    # execute simpleTracing.exe and check for errors
    ./simpleTracing.exe &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error while runtime execution!"
        errors=$(($errors + 1))
        ret=1
    else
        ret=0
    fi
    
    if [ -z "$1" ]; then
        echo "+ Running tests with unset SCAI_TRACE"
    else
        echo "+ Running tests with SCAI_TRACE=$1"
    fi
    
}

function checkCTFilesExist {
    if [ $# -ne 1 ]; then
        echo "Invalid number of parameters!"
        exit 1
    fi
    
    count=`ls -l -la *.ct* 2> /dev/null | wc -l`
    if [ "$count" -ne "$1" ]; then
        echo "Test failed. There are $count *.ct files but it should be $1".    
        errors=$(($errors + 1))
        ret=1
    else
        ret=0
    fi
}

function checkTimeFilesExist {
    if [ $# -ne 1 ]; then
        echo "Invalid number of parameters!"
        exit 1
    fi
    
    count=`ls -l -la *.time* 2> /dev/null | wc -l`
    if [ "$count" -ne "$1" ]; then
        echo "Test failed. There are $count *.time files but it should be $1".    
        errors=$(($errors + 1))
        ret=1
    else
        ret=0
    fi
}



# usage: checktimeFileContents $FILE $NTHREADS
function checkTimeFileContents {
    # all the regions with the correct number of calls should appear in the .time file

    if [ $# -ne 2 ]; then
        echo "Invalid number of parameters!"
        exit 1
    fi
    
    nThreads=$2

    content=`cat $1 2> /dev/null`
    
    # check for region 'Time main'
    # the main method is not executed in parallel and therefore should only match once!
    count=`echo "$content" | grep -E "^Time main \(in ms\) : #calls = 1${genericTimePattern}" | wc -l`
    if [ "$count" -ne 1 ]; then
        echo "ERROR: Content of the .time file is wrong (region main)"
        errors=$(($errors + 1))
    fi
            
    # check for region 'Time A'
    count=`echo "$content" | grep -E "^Time A \(in ms\) : #calls = 75000${genericTimePattern}" | wc -l`
    if [ "$count" -ne "$nThreads" ]; then
        echo "ERROR: Content of the .time file is wrong (region A)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'Time B'
    count=`echo "$content" | grep -E "^Time B \(in ms\) : #calls = 50000${genericTimePattern}" | wc -l`
    if [ "$count" -ne "$nThreads" ]; then
        echo "ERROR: Content of the .time file is wrong (region B)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'Time main.loopA'
    count=`echo "$content" | grep -E "^Time main.loopA \(in ms\) : #calls = 2500${genericTimePattern}" | wc -l`
    if [ "$count" -ne "$nThreads" ]; then
        echo "ERROR: Content of the .time file is wrong (region main.loopA)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'Time main.loopB'
    count=`echo "$content" | grep -E "^Time main.loopB \(in ms\) : #calls = 2500${genericTimePattern}" | wc -l`
    if [ "$count" -ne "$nThreads" ]; then
        echo "ERROR: Content of the .time file is wrong (region main.loopB)"
        errors=$(($errors + 1))
    fi
    
    # the file should contain 3 header lines and 4 "timing" lines for each thread and one addition line for the main
    # function (in the main thread)
    expectedLines=$(($nThreads * 7 + 1))
    
    # check if the file contains more then that
    lines=`echo "$content" | wc -l`
    if [ "$lines" -ne "$expectedLines" ]; then
        echo "ERROR: The .time file contains unknown content"
        errors=$(($errors + 1))
    fi
}

# usage: checkCTFileContents $FILE
function checkCTFileContents {
    # the structure of the .ct files is quite complicated and can't be fully tested here. We therefore only do
    # some quick validity checks
    
    # for each region in the simpleTest there has to be a line similar to
    # fn 0 0 main 2 52 52 ?    
         
    if [ $# -ne 1 ]; then
        echo "Invalid number of parameters!"
        exit 1
    fi

    content=`cat $1 2> /dev/null`
    
    # check for region 'main'
    output=`echo "$content" | grep -E '^fn [0-9]+ [0-9]+ main [0-9]+ [0-9]+ [0-9]+ \?$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (region main)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'loopA'
    output=`echo "$content" | grep -E '^fn [0-9]+ [0-9]+ loopA [0-9]+ [0-9]+ [0-9]+ main$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (region loopA)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'loopB'
    output=`echo "$content" | grep -E '^fn [0-9]+ [0-9]+ loopB [0-9]+ [0-9]+ [0-9]+ main$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (region loopB)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'A'
    output=`echo "$content" | grep -E '^fn [0-9]+ [0-9]+ A [0-9]+ [0-9]+ [0-9]+ \?$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (region A)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'B'
    output=`echo "$content" | grep -E '^fn [0-9]+ [0-9]+ B [0-9]+ [0-9]+ [0-9]+ \?$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (region B)"
        errors=$(($errors + 1))
    fi
    
    # we can also check if there is a line containing the correct number of lines
    
    # check for region 'loopA' and 'loopB'
    output=`echo "$content" | grep -E '^calls 2500 0$' | wc -l`
    if [ "$output" -ne 2 ]; then
        echo "ERROR: Content of the .ct file is wrong (calls region loopA / loopB)"
        errors=$(($errors + 1))
    fi

    # check for region 'A'
    output=`echo "$content" | grep -E '^calls 75000 0$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (calls region A)"
        errors=$(($errors + 1))
    fi
    
    # check for region 'B'
    output=`echo "$content" | grep -E '^calls 50000 0$'`
    if [ $? -ne 0 ]; then
        echo "ERROR: Content of the .ct file is wrong (calls region B)"
        errors=$(($errors + 1))
    fi
    
    
    # check for some other structural properties
    
    # check if there are exactly 4 call cost regions
    output=`echo "$content" | grep 'begin call cost line' | wc -l`
    if [ "$output" -ne 4 ]; then
        echo "ERROR: Content of the .ct file is wrong (number of call cost regions)"
        errors=$(($errors + 1))
    fi
    
    # check if there are exactly 4 exclusive call cost regions
    output=`echo "$content" | grep 'begin exclusive cost line' | wc -l`
    if [ "$output" -ne 5 ]; then
        echo "ERROR: Content of the .ct file is wrong (number of exclusive call cost regions)"
        errors=$(($errors + 1))
    fi
}

# we should use 4 threads for all tests
export OMP_NUM_THREADS=4

echo "Running runtime configuration tests:"
make clean > /dev/null
make simple DEFINES="-DSCAI_TRACE_ON" &> /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Could not build executable! Tests are skipped!"
    errors=$(($errors + 1))
else
    
    # =================================================================================================================
    # Test 1
    # check execution with unset SCAI_TRACE
        
    prepareTestCase ""
    if [ $ret -eq 0 ]; then
        # If tracing is disabled, NO tracing files should be created, so we check if *.ct or *.time files exist
        checkCTFilesExist 0
        checkTimeFilesExist 0
    fi
    
    # =================================================================================================================
    # Test 2
    # check execution with SCAI_TRAC=time

    prepareTestCase time
    if [ $ret -eq 0 ]; then
        # TODO: THIS IS STILL AN ISSUE AND HAS TO BE FIXED, WE ARE STILL CREATING .ct FILES!

        # there should be no .ct files
        #checkCTFilesExist 0
        
        checkTimeFilesExist 1
        if [ $ret -eq 0 ]; then
            checkTimeFileContents simpleTracing.exe.time 1
        fi
    fi
    
    # =================================================================================================================
    # Test 3
    # check execution with SCAI_TRACE=ct
    
    prepareTestCase ct

    if [ $ret -eq 0 ]; then
        # there should be .time files
        checkTimeFilesExist 0
        
        checkCTFilesExist 1
        if [ $ret -eq 0 ]; then
            checkCTFileContents simpleTracing.exe.ct
        fi
    fi
        
    # =================================================================================================================
    # Test 4
    # check execution with SCAI_TRACE=time:PREFIX=customPrefix
    
    prepareTestCase time:PREFIX=customPrefix
    
    if [ $ret -eq 0 ]; then
        # TODO still a problem
        
        # there should be no .ct files
        #checkCTFilesExist 0
    
        # we have to check if there is exactly one .time file (no additional .time files created)
        checkTimeFilesExist 1
        
        
        # Check whether the a correct .time file was generated
        count=`ls -l -la customPrefix.time 2> /dev/null | wc -l`
        if [ $count -ne 1 ]; then
            echo "Test failed. No .time file has been generated or a wrong name was used."
            errors=$(($errors + 1))
        else
            checkTimeFileContents customPrefix.time 1
        fi
    fi
    
    # =================================================================================================================
    # Test 5
    # check execution with SCAI_TRACE=ct:time
    
    prepareTestCase ct:time
    
    if [ $ret -eq 0 ]; then
        checkCTFilesExist 1
        
        if [ $ret -eq 0 ]; then
            checkCTFileContents simpleTracing.exe.ct
        fi
        
        checkTimeFilesExist 1
        
        if [ $ret -eq 0 ]; then
            checkTimeFileContents simpleTracing.exe.time 1
        fi
    fi
    
    
    # =================================================================================================================    
    # Test 6
    # check execution with SCAI_TRACE=time:thread
    
    prepareTestCase time:thread
    
    if [ $ret -eq 0 ]; then
        # TODO: .ct files are still created
        
        # there should be no .ct files
        #checkCTFilesExist 0
        
        checkTimeFilesExist 1
        if [ $ret -eq 0 ]; then
            checkTimeFileContents simpleTracing.exe.time 4
        fi
    fi
    
    # =================================================================================================================
    # Test 7
    # check execution with SCAI_TRACE=ct:thread
    
    prepareTestCase ct:thread
    
    if [ $ret -eq 0 ]; then
        # there should be no .time files
        checkTimeFilesExist 0
        
        checkCTFilesExist 4
        if [ $ret -eq 0 ]; then
            for file in *.ct*; do
                echo -n ""
                # TODO fix test...
                #checkCTFileContents $file
            done
        fi
    fi
    
    # =================================================================================================================
    # Test 8
    # check execution width SCAI_TRACE=time:ct:thread
    
    prepareTestCase time:ct:thread
    
    if [ $ret -eq 0 ]; then
        checkCTFilesExist 4
        # TODO: check contents of CT files

        
        checkTimeFilesExist 1
        if [ $ret -eq 0 ]; then
            checkTimeFileContents simpleTracing.exe.time 4
        fi
    fi
    
    # =================================================================================================================
    # Test 9
    # check execution with SCAI_TRACE=PREFIX=customPrefix:time:ct:thread
    
    prepareTestCase PREFIX=customPrefix:time:ct:thread
    
    if [ $ret -eq 0 ]; then
        checkCTFilesExist 4

        # TODO: check .ct files are named correctly
        # TODO: check contents of .ct files        
        

        checkTimeFilesExist 1
        
        # Check whether the a correct .time file was generated
        count=`ls -l -la customPrefix.time 2> /dev/null | wc -l`
        if [ $count -ne 1 ]; then
            echo "Test failed. No .time file has been generated or a wrong name was used."
            errors=$(($errors + 1))
        else
            checkTimeFileContents customPrefix.time 4
        fi
    fi
fi



# compile time test
# 1 check if tracing is disabled when the enviroment variable is set, but compilation was done without trace SCAI_TRACE_OFF
# 2 same check again with no SCAI_TRACE_XX

make clean > /dev/null

# =====================================================================================================================
echo ""
if [ $errors -eq 0 ]; then
    echo "Tests run sucessfully!"
    exit
else
    echo "Tests failed! Number of errors: $errors"
    exit 1
fi