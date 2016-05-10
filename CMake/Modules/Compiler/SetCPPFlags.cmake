###
 # @file SetCPPFlags.cmake
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief Set additional flags for CXX compiler and linker
 # @author Thomas Brandes
 # @date 17.07.2015
###

#### compiler dependent flag definition ####
#
#  Input variables: CMAKE_COMPILER_IS_GNUCXX, CMALE_CXX_COMPILER_ID
#
#  Output variables: 
#
#      SCAI_WARNING_FLAGS       : flags will be used only for compilation, suppress/enable warnings
#      SCAI_CXX_FLAGS           : flags used for compilation/linking
#      SCAI_CXX_FLAGS_DEBUG     : flags used for compilation/linking, only in Debug mode
#      SCAI_CXX_FLAGS_RELASE    : flags used for compilation/linking, only in Release mode
#      SCAI_CODE_COVERAGE_FLAGS : flags used for compilation/linking, only if code coverage will be enabled
#
#  All variables are only set here and not added for CMAKE variables.
#
#  Note: most common flags will be set by CMAKE itself, like -g for Debug or -O3 for Release, so these
#        flags do not have to be defined here.

set ( SCAI_LINKER_FLAGS "" )

# GNU C++ compiler

if ( CMAKE_COMPILER_IS_GNUCXX )

    set ( SCAI_WARNING_FLAGS       "-Wextra -Wall -Werror" )
    set ( SCAI_CXX_FLAGS           "" )
    set ( SCAI_CXX_FLAGS_DEBUG     "" )
    set ( SCAI_CXX_FLAGS_RELEASE   "-ffast-math -msse4a " )
    set ( SCAI_CODE_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )

    ###  Code coverage with gcov/lcov

endif ( CMAKE_COMPILER_IS_GNUCXX )


# INTEL compiler

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

    set ( MIC_NO_OFFLOAD_FLAG "-no-offload" )

    if    ( IntelCXX_COMPILER_VERSION VERSION_GREATER 14 )
        set ( MIC_NO_OFFLOAD_FLAG "-qno-offload" )
    endif ( IntelCXX_COMPILER_VERSION VERSION_GREATER 14 )

    # -fPIC should always be enabled so static libraries can be linked with shared libraries

    set ( SCAI_CXX_FLAGS "-fPIC -shared-intel " ) 

    # -wd1478 : supprress warning deprecated auto_ptr
    # not set: -Werror-all (all warnings will be errors)

    set ( SCAI_WARNING_FLAGS "-w2 -Wall -Werror-all -Wcheck -wd1478  -diag-disable 654" )
    
    # -ipo for interprocedural analysis, might be added for RELEASE but increases compile/link time dramatically
    # -xHost optimizes for the processor on which the code is compiled, not recommended for cross compilation
    #  or HPC clusters where compile node has different processor than compute nodes

    set ( SCAI_CXX_FLAGS_DEBUG "" )

    set ( SCAI_CXX_FLAGS_RELEASE "-no-prec-div " )

    #  Intel compiler requires following flags to instrument program for code coverage

    set ( SCAI_CODE_COVERAGE_FLAGS "-prof-gen=srcpos" )

endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )


# PGI C++ compiler

if ( CMAKE_CXX_COMPILER_ID MATCHES PGI )

    # Flag BOOST_HAS_THREADS is workaround needed for PGI compiler when compiling codes
    # with Boost headers

    # -DBOOST_HAS_THREADS was neded in previous versions
    # -Mipa=libc   was used in previous compiler releases
    # --gnu absolutely required if linking GNU compatible libraries, PGI has other name mangeling

    set ( SCAI_CXX_FLAGS         "-fPIC -Kieee --gnu" ) 
    set ( SCAI_CXX_FLAGS_DEBUG   "" )
    set ( SCAI_CXX_FLAGS_RELEASE "-fast " )

    # Disable warning 1097 to avoid warnings from openmpi headers with gcc specific attributes

    set ( SCAI_WARNING_FLAGS "--display_error_number --diag_suppress1097 " )

endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )
