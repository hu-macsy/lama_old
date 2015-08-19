###
 # @file LAMASummary.cmake
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
 # @brief LAMA Summary for build configuration
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

### Summary ###
message ( STATUS "" )
#message ( STATUS "        _    _")
#message ( STATUS "       ( \\__//)")
#message ( STATUS "       .'     )")
#message ( STATUS "    __/b d  .  )")
#message ( STATUS "   (_Y_`,     .)")
#message ( STATUS "    `--'-,-'  )")
#message ( STATUS "         (.  )")
#message ( STATUS "         (   )")
#message ( STATUS "        (   )")
#message ( STATUS "       ( . )         .---.")
#message ( STATUS "      (    )        (     )")
#message ( STATUS "      (   . )      (  .    )")
#message ( STATUS "      (      )    (      .  ),")
#message ( STATUS "      ( .     `\"'`  .       `)\\")
#message ( STATUS "       (      .              .)\\")
#message ( STATUS "       ((  .      .   (   .   )\\\\")
#message ( STATUS "       ((       .    (        ) \\\\")
#message ( STATUS "        ((     )     _( .   . )  \\\\")
#message ( STATUS "        ( ( .   )\"'\"`(.(     )   ( ;")
#message ( STATUS "        ( (    )      ( ( . )     \\'")
#message ( STATUS "         |~(  )        |~(  )")
#message ( STATUS "         | ||~|        | ||~|")
#message ( STATUS "    jgs  | || |        | || |  ")  
#message ( STATUS "        _| || |       _| || |")
#message ( STATUS "       /___(| |      /___(| |")
#message ( STATUS "          /___(         /___(")
message ( STATUS "" )
message ( STATUS "Summary of LAMA Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

lama_status_message ( HEADLINE "Compiler:" )
# C++ Compiler
lama_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )

# C Compiler
lama_summary_message ( "FOUND"
                       "CMAKE_C_COMPILER"
                       "C Compiler"
                       "${CMAKE_C_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CC_COMPILER_VERSION}" )

message ( STATUS "" )

lama_summary_message ( "USE"
                       "LAMA_USE_OPENMP"
                       "  OpenMP usage"
                       "" )
message ( STATUS "       OpenMP schedule set to \"${LAMA_OMP_SCHEDULE}\"" )

# LAMA (core)
message ( STATUS "" )
lama_status_message ( HEADLINE "LIBRARIES:" )

if    ( LAMA_BLAS_FOUND AND Boost_INCLUDE_DIR )
    set( REQUIRED_FOUND TRUE )
    if ( LAMA_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( LAMA_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
else  ( LAMA_BLAS_FOUND AND Boost_INCLUDE_DIR )
    set( REQUIRED_FOUND FALSE )
endif ( LAMA_BLAS_FOUND AND Boost_INCLUDE_DIR ) 

lama_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "LAMA (core)"
                       " " )
   # BLAS
    lama_summary_message ( "FOUND"
                           "LAMA_BLAS_FOUND"
                           "BLAS"
                           "(${LAMA_BLAS_NAME}) with libraries: ${LAMA_BLAS_LIBRARIES}" )
    if    ( LAMA_BLAS_NAME MATCHES "BLAS" )
    lama_summary_message ( "FOUND"
                           "LAPACK_FOUND"
                           "LAPACK"
                           "" )
    endif ( LAMA_BLAS_NAME MATCHES "BLAS" )
    
    # Boost
    lama_summary_message ( "FOUND"
                           "Boost_INCLUDE_DIR"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${Boost_INCLUDE_DIR}" )

# LAMA MPI
message ( STATUS "" )
lama_summary_message ( "USE"
                       "LAMA_USE_MPI"
                       "LAMA MPI"
                       "" )

    # MPI
    lama_summary_message ( "FOUND"
                           "MPI_FOUND"
                           "MPI"
                           "at ${MPI_INCLUDE_PATH}" )

# LAMA CUDA
message ( STATUS "" )
lama_summary_message ( "USE"
                       "LAMA_USE_CUDA"
                       "LAMA CUDA"
                       "" )

    # CUDA
    lama_summary_message ( "FOUND"
                           "CUDA_FOUND"
                           "CUDA"
                           "${CUDA_VERSION} at ${CUDA_INCLUDE_DIRS}" )
                           
    # CUDA Compute Capability
    lama_summary_message ( "FOUND"
                           "CUDA_HAVE_GPU"
                           "Compute Capability"
                           "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
message ( STATUS "" )
lama_summary_message ( "USE"
                       "LAMA_USE_MIC"
                       "LAMA MIC"
                       "" )

# LAMA TEST
message ( STATUS "" )
lama_summary_message ( "USE"
                       "LAMA_USE_CUDA"
                       "LAMA TEST"
                       "" )

    # Boost Test-Framework
    lama_summary_message ( "FOUND"
                           "Boost_UNIT_TEST_FRAMEWORK_FOUND"
                           "Boost Unit Test"
                           "" )
                           
    # Boost Regex
    lama_summary_message ( "FOUND"
                           "Boost_REGEX_FOUND"
                           "Boost Regex"
                           "" )

message ( STATUS "" )
lama_status_message ( HEADLINE "DOCUMENTATION:" )


# DOXYGEN
lama_summary_message ( "FOUND"
                       "DOXYGEN_FOUND"
                       "DOXYGEN "
                       "'make doc' to build system documentation" )
                       
# DOXYGEN
lama_summary_message ( "FOUND"
                       "SPHINX_FOUND"
                       "SPHINX"
                       "'make userdoc' to build user documentation" )

## Metis
if    ( METIS_FOUND )
    lama_status_message ( INFO "[FOUND]" "    Metis at ${METIS_INCLUDE_DIR}" )
else  ( METIS_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "Metis" )
endif ( METIS_FOUND )

# ParMetis
if    ( PARMETIS_FOUND )
    lama_status_message ( INFO "[FOUND]" "    ParMetis at ${PARMETIS_INCLUDE_DIR}" )
else  ( PARMETIS_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "ParMetis not found" )
endif ( PARMETIS_FOUND ) 

message ( STATUS "" )

lama_status_message ( HEADLINE "INFO:" )
message ( STATUS "LAMA Version : ${LAMA_VERSION} ${LAMA_VERSION_NAME}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "Library Type : ${LAMA_LIBRARY_TYPE}" )
message ( STATUS "LOG Level    : ${SCAI_LOG_LEVEL}" )
message ( STATUS "ASSERT Level : ${SCAI_ASSERT_LEVEL}" )
message ( STATUS "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
message ( STATUS "" )


# Check if all required packages are found
# LAMA (core)

if    ( NOT ( LAMA_BLAS_FOUND AND Boost_INCLUDE_DIR  ) OR ( (LAMA_BLAS_NAME MATCHES "BLAS") AND NOT LAPACK_FOUND )  )
    message( FATAL_ERROR "Configuration for LAMA (core) incomplete!")
endif ( NOT ( LAMA_BLAS_FOUND AND Boost_INCLUDE_DIR  ) OR ( (LAMA_BLAS_NAME MATCHES "BLAS") AND NOT LAPACK_FOUND ) )

# LAMA MPI
if    ( LAMA_USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "Configuration for LAMA MPI incomplete!")
endif ( LAMA_USE_MPI AND NOT MPI_FOUND )

# LAMA MPI
if    ( LAMA_USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "Build of LAMA MPI enabled, but configuration is incomplete!")
endif ( LAMA_USE_MPI AND NOT MPI_FOUND )

# LAMA Cuda
if    ( LAMA_USE_CUDA AND NOT CUDA_FOUND )
    message( FATAL_ERROR "Build of LAMA Cuda enabled, but configuration is incomplete!")
endif ( LAMA_USE_CUDA AND NOT CUDA_FOUND )

# LAMA Test
if    ( LAMA_BUILD_TEST AND NOT ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND ) )
    message( FATAL_ERROR "Build of LAMA Test enabled, but configuration is incomplete!")
endif ( LAMA_BUILD_TEST AND NOT ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND ) )
