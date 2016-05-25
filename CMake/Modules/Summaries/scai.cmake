###
 # @file Summaries/scai.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief scai Summary for build configuration
 # @author Jan Ecker
 # @date 25.04.2013
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )

include ( Summaries/Modules/Compiler )

#lama all core
heading ( "Required core:" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )
    set ( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )

heading2 ( "External Libraries" "REQUIRED_FOUND" )

    # pthreads
    found_message ( "pThreads" "SCAI_THREAD_LIBRARIES" "REQUIRED" "Version ${SCAI_THREAD_VERSION}" )
    # boost
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

    include ( Summaries/Modules/BLAS )

heading ( "Optional External Libraries:" )
include ( Summaries/Modules/Accelerator )
include ( Summaries/Modules/Distributed )
include ( Summaries/Modules/Graphpartitioning )

heading ( "Optional components:" "" )
heading3 ( "Java:" "JAVA_FOUND" )
    found_message ( "Java Exexutable" "JAVA_FOUND" "OPTIONAL" "with ${Java_JAVAC_EXECUTABLE}" )

include ( Summaries/Modules/Build )

include ( Summaries/Modules/Configuration )
