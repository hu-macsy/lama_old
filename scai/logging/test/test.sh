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

genericPatternSimple=[0-9]{4}-[0-9]{2}-[0-9]{2},\ [0-9]{2}:[0-9]{2}:[0-9]{2}" Test @ <unk_thread> \( main -> simpleLogging.cpp::"[0-9]{1,2}" \)"
genericPatternComplex1=[0-9]{4}-[0-9]{2}-[0-9]{2},\ [0-9]{2}:[0-9]{2}:[0-9]{2}
genericPatternComplex2="@ <unk_thread> \( main -> complexLogging.cpp::"[0-9]{1,2}" \)"
errors=0

# =====================================================================================================================
# Runtime configurations tests
#
# The simple executable is build with full logging (level TRACE). All different log levels are set using export and the
# output is checked.
# =====================================================================================================================
echo "Running runtime configuration tests:"
make clean > /dev/null
make simple DEFINES="-DSCAI_LOG_LEVEL_TRACE" > /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Could not build executable! Tests are skipped!"
    errors=$(($errors + 1))
else
    for level in "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "FATAL" "OFF"; do
	export SCAI_LOG=$level
	output=$( ./simpleLogging.exe | tr -d '\n' ) 

	pattern=^
	case $level in
	"TRACE")
	    pattern+=${genericPatternSimple}" TRACE trace message"
	    ;&
	"DEBUG")
	    pattern+=${genericPatternSimple}" DEBUG debug message"
	    ;&
	"INFO")
	    pattern+=${genericPatternSimple}" INFO info message"
	    ;&
	"WARN")
	    pattern+=${genericPatternSimple}" WARN warn message"
	    ;&
	"ERROR")
	    pattern+=${genericPatternSimple}" ERROR error message"
	    ;&
	"FATAL")
	    pattern+=${genericPatternSimple}" FATAL fatal message"
	    ;&
	esac
	pattern+=$

	if ! [[ "$output" =~ $pattern ]]; then
	    echo "ERROR: Output of the log level $level is wrong."
	    errors=$(($errors + 1))
	fi
    done
    echo "done"
fi

# =====================================================================================================================
# Compile-time configuration tests
#
# The exported debug level is set to full logging (level TRACE). The simple executable is build using all available log
# levels and the ouput is checked.
# =====================================================================================================================
echo ""
echo "Running compile-time configuration tests:"
export SCAI_LOG=TRACE

for level in "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "FATAL" "OFF"; do
    make clean > /dev/null
    make simple DEFINES="-DSCAI_LOG_LEVEL_${level}" > /dev/null
    if [ $? -ne 0 ]; then
	echo "ERROR: Could not build executable using log level $level"
	errors=$(($errors + 1))
    else
	output=$( ./simpleLogging.exe | tr -d '\n' ) 

	pattern=^
	case $level in
	"TRACE")
	    pattern+=${genericPatternSimple}" TRACE trace message"
	    ;&
	"DEBUG")
	    pattern+=${genericPatternSimple}" DEBUG debug message"
	    ;&
	"INFO")
	    pattern+=${genericPatternSimple}" INFO info message"
	    ;&
	"WARN")
	    pattern+=${genericPatternSimple}" WARN warn message"
	    ;&
	"ERROR")
	    pattern+=${genericPatternSimple}" ERROR error message"
	    ;&
	"FATAL")
	    pattern+=${genericPatternSimple}" FATAL fatal message"
	    ;&
	esac
	pattern+=$

	if ! [[ "$output" =~ $pattern ]]; then
	    echo "ERROR: Output of the log level $level is wrong."
	    errors=$(($errors + 1))
	fi
    fi
done
echo "done"
# =====================================================================================================================
# Hierarchical configuration tests
#
# The complex executable is build with full logging (level TRACE). Different log levels for the different loggers in
# the hierarchy are set and the output is checked.
# =====================================================================================================================
echo ""
echo "Hierarchical configuration tests:"
make clean > /dev/null
make complex DEFINES="-DSCAI_LOG_LEVEL_TRACE" > /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Could not build executable! Tests are skipped!"
    errors=$(($errors + 1))
else
    export SCAI_LOG=loggerConfig.cfg

    # Check change of root debug level by reducing default from INFO to DEBUG
    echo "<root>=DEBUG" > loggerConfig.cfg
    
    # Check simple change of log level on second level
    echo "Class1.method1=TRACE" >> loggerConfig.cfg

    # Check use of three levels with different log levels
    echo "Class1.method2=OFF" >> loggerConfig.cfg
    echo "Class1.method2.region1=DEBUG" >> loggerConfig.cfg

    # Check inheritance of log levels when level is changed on the first level
    echo "Class2=WARN" >> loggerConfig.cfg
    echo "Class2.method2=TRACE" >> loggerConfig.cfg

    output=$( ./complexLogging.exe | tr -d '\n' )

    pattern=^
    pattern+=$genericPatternComplex1" Class1 "$genericPatternComplex2" DEBUG message class1"
    pattern+=$genericPatternComplex1" Class1.method1 "$genericPatternComplex2" TRACE message class1 method1"
    pattern+=$genericPatternComplex1" Class1.method2.region1 "$genericPatternComplex2" INFO message class1 method2 region1"
    pattern+=$genericPatternComplex1" Class2 "$genericPatternComplex2" WARN message class2"
    pattern+=$genericPatternComplex1" Class2.method2 "$genericPatternComplex2" TRACE message class2 method2"
    pattern+=$

    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Wrong output!"
	errors=$(($errors + 1))
    fi

    rm loggerConfig.cfg
    echo "done"
fi

# =====================================================================================================================
# Format string configuration tests
# 
# The simple executable is build with log level FATAL only. The output format is changed using a config file and the
# output is checked.
# =====================================================================================================================
echo ""
echo "Format string configuration tests:"
make clean > /dev/null
make simple DEFINES="-DSCAI_LOG_LEVEL_FATAL" > /dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: Could not build executable! Tests are skipped!"
    errors=$(($errors + 1))
else
    export SCAI_LOG=loggerConfig.cfg

    # check all available variables
    # check format string: #date
    echo 'format="#date"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^[0-9]{4}-[0-9]{2}-[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #date"
	errors=$(($errors + 1))
    fi

    # check format string: #time
    echo 'format="#time"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^[0-9]{2}:[0-9]{2}:[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #time"
	errors=$(($errors + 1))
    fi

    # check format string: #name
    echo 'format="#name"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"Test"$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #name"
	errors=$(($errors + 1))
    fi

    # check format string: #func
    echo 'format="#func"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"main"$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #func"
	errors=$(($errors + 1))
    fi

    # check format string: #file
    echo 'format="#file"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"simpleLogging.cpp"$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #file"
	errors=$(($errors + 1))
    fi

    # check format string: #line
    echo 'format="#line"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^[0-9]{1,2}$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #line"
	errors=$(($errors + 1))
    fi

    # check format string: #level
    echo 'format="#level"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^FATAL$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #level"
	errors=$(($errors + 1))
    fi

    # check format string: #msg
    echo 'format="#msg"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"fatal message"$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using format string #msg"
	errors=$(($errors + 1))
    fi


    # check some combinations of variables
    # check combined format string: #msg #level #line
    echo 'format="#msg #level #line"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"fatal message FATAL "[0-9]{1,2}$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using combined format string #msg #level #line"
	errors=$(($errors + 1))
    fi

    # check combined format string: #file #msg #date
    echo 'format="#file #msg #date"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"simpleLogging.cpp fatal message "[0-9]{4}-[0-9]{2}-[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using combined format string #file #msg #date"
	errors=$(($errors + 1))
    fi

    # check format string with additional text in it
    echo 'format="error message: #msg"' > loggerConfig.cfg
    output=$( ./simpleLogging.exe | tr -d '\n' )
    pattern=^"error message: fatal message"$
    if ! [[ "$output" =~ $pattern ]]; then
	echo "ERROR: Output wrong when using additional text in format string"
	errors=$(($errors + 1))
    fi

    rm loggerConfig.cfg
    echo "done"
fi

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