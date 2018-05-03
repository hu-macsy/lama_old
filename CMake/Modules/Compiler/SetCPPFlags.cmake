###
 # @file SetCPPFlags.cmake
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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

if ( CMAKE_CXX_COMPILER_ID MATCHES GNU )

    set ( SCAI_WARNING_FLAGS       "-Wextra -Wall -Werror -Wno-deprecated-declarations -Wno-unknown-pragmas" )
    set ( SCAI_CXX_FLAGS           "" )
    set ( SCAI_CXX_FLAGS_DEBUG     "" )
    set ( SCAI_CXX_FLAGS_RELEASE   "-ffast-math -msse4a " )
    set ( SCAI_CODE_COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
    set ( SCAI_STATIC_FLAGS        "-fPIC" )

    ###  Code coverage with gcov/lcov

endif ( CMAKE_CXX_COMPILER_ID MATCHES GNU )


# INTEL compiler

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

    # -fPIC should always be enabled so static libraries can be linked with shared libraries

    set ( SCAI_STATIC_FLAGS "-fPIC -shared-intel " ) 

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

# Clang Compiler (llvm or AppleClang)

if ( CMAKE_CXX_COMPILER_ID MATCHES Clang )

    # "-Weverything" # wow this creates so much warnings

    set ( SCAI_WARNING_FLAGS "-Wall -Werror -Wno-deprecated-declarations -Wno-unknown-pragmas" ) 

    if ( CXX_SUPPORTS_C11 )
        set ( SCAI_CXX_FLAGS           "-stdlib=libc++" )
    endif ( CXX_SUPPORTS_C11 )

    set ( SCAI_CXX_FLAGS_DEBUG     "" )
    set ( SCAI_CXX_FLAGS_RELEASE   "-ffast-math" )
    #set ( SCAI_CODE_COVERAGE_FLAGS "-fsanitize-coverage=???" )

endif ( CMAKE_CXX_COMPILER_ID MATCHES Clang )


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

    set ( SCAI_WARNING_FLAGS "--display_error_number --diag_suppress1097 " ) #-Werror 

endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )
