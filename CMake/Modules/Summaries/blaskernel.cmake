###
 # @file Summaries/blaskernel.cmake
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
 # @brief SCAI blaskernel Configuration Summary
 # @author Eric Schricker
 # @date 11.01.2016
 # @since 2.0.0
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "=========================================" )
message ( STATUS "Summary of SCAI blaskernel Configuration:" )
message ( STATUS "=========================================" )

include ( Summaries/Modules/Compiler )

# blaskernel (core)
set ( REQUIRED_FOUND FALSE )
if    ( SCAI_BLAS_FOUND )
    set( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_BLAS_FOUND )

heading2 ( "Required core" "REQUIRED_FOUND" )

    # BLAS (Lapack)
    found_message ( "BLAS" "SCAI_BLAS_FOUND" "REQUIRED" "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )
    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
        message ( STATUS "                                 ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        found_message ( "Lapack" "LAPACK_FOUND" "REQUIRED" "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND )

heading2 ( "Optional components" "" )
include ( Summaries/Modules/Accelerator )

heading3 ( "Internal Libraries" "REQUIRED_FOUND" )
    found_message ( "SCAI common"    "SCAI_COMMON_FOUND"    "REQUIRED" "Version ${SCAI_COMMON_VERSION}"    )
    found_message ( "SCAI logging"   "SCAI_LOGGING_FOUND"   "REQUIRED" "Version ${SCAI_LOGGING_VERSION}"   )
    found_message ( "SCAI tracing"   "SCAI_TRACING_FOUND"   "REQUIRED" "Version ${SCAI_TRACING_VERSION}"   )
    found_message ( "SCAI tasking"   "SCAI_TASKING_FOUND"   "REQUIRED" "Version ${SCAI_TASKING_VERSION}"   )
    found_message ( "SCAI kregistry" "SCAI_KREGISTRY_FOUND" "REQUIRED" "Version ${SCAI_KREGISTRY_VERSION}" )

include ( Summaries/Modules/Build )
                       
include ( Summaries/Modules/Configuration )
