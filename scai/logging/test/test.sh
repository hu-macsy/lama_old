#!/bin/bash
###
 # @file scai/logging/test/test.sh
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
 # @brief Tests for SCAI logging
 # @author Jan Ecker
 # @date 03.09.2015
###

if ((BASH_VERSINFO[0] < 4))
then
	echo "For testing logging you need bash version 4 or newer"
	exit
fi

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BINDIR=$SCRIPTDIR

genericPatternSimple=[0-9]{4}-[0-9]{2}-[0-9]{2},\ [0-9]{2}:[0-9]{2}:[0-9]{2}" Test @ thread_"[0-9]{1,2}" \( main -> simpleLogging.cpp::"[0-9]{1,2}" \)"
genericPatternComplex1=[0-9]{4}-[0-9]{2}-[0-9]{2},\ [0-9]{2}:[0-9]{2}:[0-9]{2}
genericPatternComplex2="@ thread_"[0-9]{1,2}" \( main -> complexLogging.cpp::"[0-9]{1,2}" \)"
errors=0

# =====================================================================================================================
# Runtime configurations tests
#
# The simple executable is build with full logging (level TRACE). All different log levels are set using export and the
# output is checked.
# =====================================================================================================================
echo "Running runtime configuration tests:"

for level in "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "FATAL" "OFF"; do
    export SCAI_LOG=$level
    output=$( $BINDIR/simpleLoggingTRACE | tr -d '\n' )

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

# =====================================================================================================================
# Compile-time configuration tests
#
# The exported debug level is set to full logging (level TRACE). The simple executable is build using all available log
# levels and the output is checked.
# =====================================================================================================================
echo ""
echo "Running compile-time configuration tests:"
export SCAI_LOG=TRACE

for level in "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "FATAL" "OFF"; do

        output=$( $BINDIR/simpleLogging${level} | tr -d '\n' )

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

# =====================================================================================================================
# Hierarchical configuration tests
#
# The complex executable is build with full logging (level TRACE). Different log levels for the different loggers in
# the hierarchy are set and the output is checked.
# =====================================================================================================================
echo ""
echo "Hierarchical configuration tests:"

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

    output=$( $BINDIR/complexLogging | tr -d '\n' )

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

# =====================================================================================================================
# Format string configuration tests
# 
# The simple executable is build with log level FATAL only. The output format is changed using a config file and the
# output is checked.
# =====================================================================================================================

echo ""
echo "Format string configuration tests:"

    export SCAI_LOG=loggerConfig.cfg

    # check all available variables
    # check format string: #date
    echo 'format="#date"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^[0-9]{4}-[0-9]{2}-[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #date"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #time
    echo 'format="#time"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^[0-9]{2}:[0-9]{2}:[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #time"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #name
    echo 'format="#name"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"Test"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #name"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #func
    echo 'format="#func"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"main"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #func"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #file
    echo 'format="#file"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"simpleLogging.cpp"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #file"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #line
    echo 'format="#line"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^[0-9]{1,2}$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #line"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #level
    echo 'format="#level"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^FATAL$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #level"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string: #msg
    echo 'format="#msg"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"fatal message"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #msg"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi
    
    # check format string: #stack
    echo 'format="#stack"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern="stack"
    # Only for GCC
    #pattern=^"    stack\[1\] : scai::logging::GenLogger::log\(char const\*, scai::logging::SourceLocation&, std::string const&\)"
    #pattern+="    stack\[2\] : scai::logging::GenLogger::fatal\(scai::logging::SourceLocation, std::string const&\)"
    #pattern+=.*"\[0x"[0-9a-f]+"\]"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #stack"
        errors=$(($errors + 1))
    fi
    
    # check format string: #<environmental variable> with existing variable
    echo 'format="#ENV_TEST_VAR"' > loggerConfig.cfg
    export ENV_TEST_VAR=scai_test
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"scai_test"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #<environmental variable> with existing variable"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi
    
    # check format string: #<environmental variable> with not existing variable
    echo 'format="#ENV_TEST_VAR"' > loggerConfig.cfg
    unset ENV_TEST_VAR
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"\\$\{ENV_TEST_VAR\}"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #<environmental variable> with not existing variable"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi
    
    # check format string: #MsG (check that variable parsing is case insensitive)
    echo 'format="#msg"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"fatal message"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using format string #MsG (variables case sensitive?)"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check some combinations of variables
    # check combined format string: #msg #level #line
    echo 'format="#msg #level #line"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"fatal message FATAL "[0-9]{1,2}$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using combined format string #msg #level #line"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check combined format string: #file #msg #date
    echo 'format="#file #msg #date"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"simpleLogging.cpp fatal message "[0-9]{4}-[0-9]{2}-[0-9]{2}$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using combined format string #file #msg #date"
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi

    # check format string with additional text in it
    echo 'format="error message: #msg"' > loggerConfig.cfg
    output=$( $BINDIR/simpleLoggingFATAL | tr -d '\n' )
    pattern=^"error message: fatal message"$
    if ! [[ "$output" =~ $pattern ]]; then
        echo "ERROR: Output wrong when using additional text in format string."
        echo "Looking for pattern: \"$pattern\" but output was \"$output\""
        errors=$(($errors + 1))
    fi
    
    

    rm loggerConfig.cfg
    echo "done"


# =====================================================================================================================
echo ""
if [ $errors -eq 0 ]; then
    echo "Tests run sucessfully!"
    exit
else
    echo "Tests failed! Number of errors: $errors"
    exit 1
fi
