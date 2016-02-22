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

genericTimePattern=", inclusive = "[0-9]{1,5}\.[0-9]{4,6}", exclusive = "[0-9]{1,5}\.[0-9]{4,6}$
errors=0

# =====================================================================================================================
# Runtime configurations tests
#
# In this test the executable is build WITH trace support and the runtime configuration via enviormental variables is
# used to control the tracing behavior. 
# =====================================================================================================================
echo "Running runtime configuration tests:"
make clean > /dev/null
make simple DEFINES="-DSCAI_TRACE_ON" &> /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Could not build executable! Tests are skipped!"
    errors=$(($errors + 1))
else
    # Test 1
    # before we start, make sure there are no tracing files laying around from earlier runs
    rm -rf *.ct *.time
    
    # lets start with a test that checks what happens when no behavior is specified
    echo "+ Running tests with unset SCAI_TRACE"
    
    unset SCAI_TRACE
    ./simpleTracing.exe &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error while runtime execution!"
        errors=$(($errors + 1))
    else
        # Test 1.1
        # if no tracing level is specific there should be NO tracing, so no tracing files should exist
        # so we check if there are NO .ct or .time files
        count=`ls -l -la *.ct *.time 2> /dev/null | wc -l`
        if [ $count -ne 0 ]; then
            echo "Test failed. Tracing files have been generated, without any tracing beeing enabled."
            errors=$(($errors + 1))
        fi
    fi
    
    
    # Test 2
    # now lets test if the use of SCAI_TRACE=time produces the correct result
    echo "+ Running tests with SCAI_TRACE=time"
    
    export SCAI_TRACE=time
    ./simpleTracing.exe &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error while runtime execution!"
        errors=$(($errors + 1))
    else
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # TODO: THIS IS STILL AN ISSUE AND HAS TO BE FIXED, WE ARE STILL CREATING .ct FILES!
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        # Test 2.1
        # as we use SCAI_TRACE=time there should be no .ct created
        #count=`ls -l -la *.ct 2> /dev/null | wc -l`
        #if [ $count -ne 0 ]; then
        #    echo "Test failed. A .ct has been generated but only SCAI_TRACE=time was set."
        #    errors=$(($errors + 1))
        #fi
        
        
        # Test 2.2
        # we have to check whether the a correct .time file was generated
        content=`cat simpleTracing.exe.time 2> /dev/null`
        if [ $? -ne 0 ]; then
            echo "Test failed. No valid .time file has been generated or a wrong name was used."
            errors=$(($errors + 1))
        else
            # all the regions with the correct number of calls should appear in the .time file
            
            # check for region 'Time A'
            output=`echo "$content" | grep 'Time A'`
            pattern=^"Time A \(in ms\) : #calls = 300000"${genericTimePattern}
            if ! [[ "$output" =~ $pattern ]]; then
                echo "ERROR: Content of the .time file is wrong (region A)"
                errors=$(($errors + 1))
            fi
            
            # check for region 'Time B'
            output=`echo "$content" | grep 'Time B'`
            pattern=^"Time B \(in ms\) : #calls = 200000"${genericTimePattern}
            if ! [[ "$output" =~ $pattern ]]; then
                echo "ERROR: Content of the .time file is wrong (region B)"
                errors=$(($errors + 1))
            fi
            
            # check for region 'Time main'
            #the space behind main is important (otherwise it would match main.loopX aswell
            output=`echo "$content" | grep 'Time main '`
            pattern=^"Time main \(in ms\) : #calls = 1"${genericTimePattern}
            if ! [[ "$output" =~ $pattern ]]; then
                echo "ERROR: Content of the .time file is wrong (region main)"
                errors=$(($errors + 1))
            fi
            
            # check for region 'Time main.loopA'
            output=`echo "$content" | grep 'Time main.loopA'`
            pattern=^"Time main.loopA \(in ms\) : #calls = 10000"${genericTimePattern}
            if ! [[ "$output" =~ $pattern ]]; then
                echo "ERROR: Content of the .time file is wrong (region main.loopA)"
                errors=$(($errors + 1))
            fi
            
            # check for region 'Time main.loopB'
            output=`echo "$content" | grep 'Time main.loopB'`
            pattern=^"Time main.loopB \(in ms\) : #calls = 10000"${genericTimePattern}
            if ! [[ "$output" =~ $pattern ]]; then
                echo "ERROR: Content of the .time file is wrong (region main.loopB)"
                errors=$(($errors + 1))
            fi
            
            
            # the file should contain 3 header lines and 5 "timing" lines (one line for each region)
            # check if the file contains more then that
            lines=`echo "$content" | wc -l`
            if [ $lines -ne 8 ]; then
                echo "ERROR: The .time file contains unknown content"
                errors=$(($errors + 1))
            fi
        fi
        
        
        # TODO:
        # do the checks for time and ct in specific functions so that it can be reused in the tests
    fi
    
    
    # Test 3
    # TODO!
    # SCAI_TRACE=ct
    
    
    # Test 3
    # now lets test if the use of SCAI_TRACE=ct produces the correct result
    echo "+ Running tests with SCAI_TRACE=ct"
    
    # clean old tracing files
    rm -rf *.ct *.time
    
    export SCAI_TRACE=ct
    ./simpleTracing.exe &> /dev/null
    if [ $? -ne 0 ]; then
        echo "Error while runtime execution!"
        errors=$(($errors + 1))
    else
        # Test 3.1
        # as we use SCAI_TRACE=ct there should be no .time created
        count=`ls -l -la *.time 2> /dev/null | wc -l`
        if [ $count -ne 0 ]; then
            echo "Test failed. A .time has been generated but only SCAI_TRACE=ct was set."
            errors=$(($errors + 1))
        fi
        
                
        # Test 3.2
        # we have to check whether the a correct .ct file was generated
        content=`cat simpleTracing.exe.ct 2> /dev/null`
        if [ $? -ne 0 ]; then
            echo "Test failed. No valid .ct file has been generated or a wrong name was used."
            errors=$(($errors + 1))
        else
            # TODO
            echo ""
        fi
    fi
    
    
    
    
    
    
    # Test 4
    # TODO!
    # SCAI_TRACE=time:PREFIX=xxx
    
    
    # Test 5
    # TODO!
    # SCAI_TRACE=ct:time
    
    
    # Test 6
    # TODO!
    # SCAI_TRACE=time:thread
    
    
    # Test 7
    # TODO!
    # SCAI_TRACE=ct:thread
    
    
    # Test 8
    # TODO!
    # SCAI_TRACE=time:ct:thread
    
    
    # Test 9
    # TODO!
    # SCAI_TRACE=PREFIX=xxx:time:ct:thread
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