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

message ( STATUS "${CMAKE_CXX_COMPILER_ID} compiler" )

#### Check for -std=c++11

include ( CheckCXXCompilerFlag )

if ( NOT DEFINED CXX_SUPPORTS_C11 )
    CHECK_CXX_COMPILER_FLAG( -std=c++11 CXX_SUPPORTS_C11 )
endif ()

#### OpenMP ####

find_package ( OpenMP )

#### Compile/Link flag for OpenMP will be set for all source files and all targets

if ( OPENMP_FOUND )
   set ( LAMA_CXX_FLAGS "${OpenMP_CXX_FLAGS}" )
else ( OPENMP_FOUND )
   set ( LAMA_CXX_FLAGS "" )
endif ( OPENMP_FOUND )

#### compiler dependent flag definition ####

# GNU
if ( CMAKE_COMPILER_IS_GNUCXX )

    if ( CXX_SUPPORTS_C11 )
        set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} -std=c++11" )
    endif ( CXX_SUPPORTS_C11 )

    set ( LAMA_LINKER_FLAGS "-Wl,--no-as-needed " )
    set ( LAMA_WARNING_FLAGS "-Wextra -Wall -Werror" ) # -pedantic -std=c++98 " ) # -march=core02

    # Supress unknown pragma warnings if OpenMP is disabled
    if ( NOT OPENMP_FOUND )
        set ( LAMA_WARNING_FLAGS "${LAMA_WARNING_FLAGS} -Wno-unknown-pragmas" )
    endif ( NOT OPENMP_FOUND )
    
    set ( LAMA_CXX_FLAGS_RELEASE "-ffast-math -msse4a " )

endif ( CMAKE_COMPILER_IS_GNUCXX )


# INTEL
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

    message ( STATUS "LAMA_CXX_FLAGS = ${LAMA_CXX_FLAGS}" )

    set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} -fPIC -shared-intel" ) 

    if ( CXX_SUPPORTS_C11 )
        set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} -std=c++11" )
    else ( CXX_SUPPORTS_C11 )
        set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} -std=c++0x" )
    endif ( CXX_SUPPORTS_C11 )

    message ( STATUS "LAMA_CXX_FLAGS = ${LAMA_CXX_FLAGS}" )
    
    # -wd1478 : supprress warning deprecated auto_ptr
    # not set: -Werror-all (all warnings will be errors)

    set ( LAMA_WARNING_FLAGS "-w2 -Wall -Wcheck -wd1478" ) # -Werror-all Warnings/Errors. No Remarks.
    
    # Supress unknown pragma warnings if OpenMP is disabled
    if ( NOT OPENMP_FOUND )
        set ( LAMA_WARNING_FLAGS "${LAMA_WARNING_FLAGS} -Wno-unknown-pragmas" )
    endif ( NOT OPENMP_FOUND )
    
    set ( LAMA_CXX_FLAGS_RELEASE "-ipo -no-prec-div -xHost " )

endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )


# PGI
if ( CMAKE_CXX_COMPILER_ID MATCHES PGI )

    set ( LAMA_CXX_FLAGS "-fPIC -Kieee -Mipa=libc -DBOOST_HAS_THREADS " ) # -std=c++0x

    # Disable warning 1097 to avoid warnings from openmpi headers with gcc specific attributes

    set ( LAMA_WARNING_FLAGS "--display_error_number --diag_suppress1097 " )
    
    set ( LAMA_CXX_FLAGS_RELEASE "-fast " )

endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )


## add variables to cache with new names so they can be modified by the user via CCMAKE

set ( ADDITIONAL_CXX_FLAGS "${LAMA_CXX_FLAGS}" CACHE STRING "additional flags for cxx compile and link" )
set ( ADDITIONAL_WARNING_FLAGS "${LAMA_WARNING_FLAGS}" CACHE STRING "compilation flags concerning warnings" )
set ( ADDITIONAL_CXX_FLAGS_RELEASE "${LAMA_CXX_FLAGS_RELEASE}" CACHE STRING "addtional cxx compiler flags for release optimizations" )
set ( ADDITIONAL_LINKER_FLAGS "${LAMA_LINKER_FLAGS}" CACHE STRING "additional linker flags" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS ADDITIONAL_WARNING_FLAGS ADDITIONAL_CXX_FLAGS_RELEASE ADDITIONAL_LINKER_FLAGS )

#### concluding all defined compiler flags to CMAKE_..._FLAGS ####

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}")
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
