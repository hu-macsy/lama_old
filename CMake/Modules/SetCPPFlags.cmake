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

# GNU
if ( CMAKE_COMPILER_IS_GNUCXX )

    #set ( LAMA_LINKER_FLAGS "-Wl,--no-as-needed " )
    set ( SCAI_WARNING_FLAGS "-Wextra -Wall -Werror" ) # -pedantic -std=c++98 " ) # -march=core02

    set ( LAMA_CXX_FLAGS_RELEASE "-ffast-math -msse4a " )

endif ( CMAKE_COMPILER_IS_GNUCXX )


# INTEL
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

    message ( STATUS "LAMA_CXX_FLAGS = ${LAMA_CXX_FLAGS}" )

    # -fPIC should always be enabled so static libraries can be linked with shared libraries

    set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} -fPIC -shared-intel " ) 

    # -wd1478 : supprress warning deprecated auto_ptr
    # not set: -Werror-all (all warnings will be errors)

    set ( SCAI_WARNING_FLAGS "-w2 -Wall -Wcheck -wd1478" ) # -Werror-all Warnings/Errors. No Remarks.
    
    # -ipo for interprocedural analysis, might be added for RELEASE but increases compile/link time dramatically
    # -xHost optimizes for the processor on which the code is compiled, not recommended for cross compilation
    #  or HPC clusters where compile node has different processor than compute nodes

    set ( LAMA_CXX_FLAGS_RELEASE "-no-prec-div " )

endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )


# PGI
if ( CMAKE_CXX_COMPILER_ID MATCHES PGI )

    # Flag BOOST_HAS_THREADS is workaround needed for PGI compiler when compiling codes
    # with Boost headers

    set ( LAMA_CXX_FLAGS "-fPIC -Kieee -Mipa=libc -DBOOST_HAS_THREADS " ) # -std=c++0x

    # Disable warning 1097 to avoid warnings from openmpi headers with gcc specific attributes

    set ( SCAI_WARNING_FLAGS "--display_error_number --diag_suppress1097 " )
    
    set ( LAMA_CXX_FLAGS_RELEASE "-fast " )

endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )

###  Code coverage with gcov/lcov
if    ( USE_CODE_COVERAGE )
    set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
    set ( LAMA_CXX_FLAGS ${LAMA_CXX_FLAGS} ${COVERAGE_FLAGS} )
endif ( USE_CODE_COVERAGE )
