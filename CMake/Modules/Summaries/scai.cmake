###
 # @file Summaries/scai.cmake
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
 # @brief scai Summary for build configuration
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

include ( Functions/scaiStatusMessage )
include ( Functions/scaiSummaryMessage )

message ( STATUS "" )
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

scai_status_message ( HEADLINE "Compiler:" "" )
message ( STATUS "" )

if    ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND TRUE )
else  ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND FALSE )
endif ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "REQUIRED"
                       "Configuration"
                       "" )

# C++ Compiler
scai_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "REQUIRED"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )

scai_summary_message ( "FOUND"
					             "CXX_SUPPORTS_C11"
                       "REQUIRED"
					             "with C++11 support"
					             "" )

if    ( NOT CXX_SUPPORTS_C11 )
    message ( STATUS "" )
    message ( STATUS "Either compiler supporting C++11 or Boost needed." )

    scai_summary_message ( "FOUND"
                           "SCAI_BOOST_INCLUDE_DIR"
                           "REQUIRED"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )

# LAMA (core)
message ( STATUS "" )
scai_status_message ( HEADLINE "External Libraries:" )
message ( STATUS "" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )
    set ( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "REQUIRED"
                       "  Required core"
                       "" )

message ( STATUS "" )

    # pthreads
    scai_summary_message ( "FOUND"
                           "SCAI_THREAD_LIBRARIES"
                           "REQUIRED"
                           "pThreads"
                           "VERSION ${SCAI_THREAD_VERSION}" )

    # boost
    scai_summary_message ( "FOUND"
                           "SCAI_BOOST_INCLUDE_DIR"
                           "REQUIRED"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

    # BLAS
    scai_summary_message ( "FOUND"
                           "SCAI_BLAS_FOUND"
                           "REQUIRED"
                           "BLAS"
                           "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )

    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
      message ( STATUS "                              ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
                           
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        scai_summary_message ( "FOUND"
                               "LAPACK_FOUND"
                               "REQUIRED"
                               "LAPACK"
                               "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )

message ( STATUS "" )
message ( STATUS "  Optional components" )
message ( STATUS "" )

# OpenMP usage

scai_summary_message ( "USE"
                       "USE_OPENMP"
                       "OPTIONAL"
                       "OpenMP"
                       ""   )

    scai_summary_message ( "FOUND"
                           "OPENMP_VERSION"
                           "OPTIONAL"
                           "OpenMP"
                           "Version ${OPENMP_VERSION}" )

    scai_summary_message ( "FOUND"
                           "OpenMP_CXX_FLAGS"
                           "OPTIONAL"
                           "compile flag"
                           "${OpenMP_CXX_FLAGS}" )

    scai_summary_message ( "FOUND"
                           "SCAI_OMP_SCHEDULE"
                           "OPTIONAL"
                           "schedule type"
                           "set to \"${SCAI_OMP_SCHEDULE}\"" )

# LAMA CUDA
set ( REQUIRED_FOUND FALSE )
if    ( CUDA_FOUND AND USE_CUDA )
  set ( REQUIRED_FOUND TRUE )
endif ( CUDA_FOUND AND USE_CUDA )

message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "CUDA"
                       "" )

    # CUDA
    scai_summary_message ( "FOUND"
                           "CUDA_FOUND"
                           "OPTIONAL"
                           "CUDA"
                           "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
                           
    # CUDA Compute Capability
    scai_summary_message ( "FOUND"
                           "CUDA_HAVE_GPU"
                           "OPTIONAL"
                           "Compute Capability"
                           "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
message ( STATUS "" )
scai_summary_message ( "USE"
                       "USE_MIC"
                       "OPTIONAL"
                       "MIC"
                       "" )

# LAMA MPI
set ( REQUIRED_FOUND FALSE )
if    ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )
  set ( REQUIRED_FOUND TRUE )
endif ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )

message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "Distributed"
                       "" )

    # MPI
    scai_summary_message ( "FOUND"
                           "MPI_FOUND"
                           "OPTIONAL"
                           "MPI"
                           "Version ${MPI_VERSION} at ${SCAI_MPI_INCLUDE_DIR}" )

    # GPI
    scai_summary_message ( "FOUND"
                           "GPI_FOUND"
                           "OPTIONAL"
                           "GPI"
                           "at ${SCAI_GPI_INCLUDE_DIR}" )

# Graph Partitioning
set ( REQUIRED_FOUND FALSE )
if    ( METIS_FOUND AND USE_GRAPHPARTITIONING )
  set ( REQUIRED_FOUND TRUE )
endif ( METIS_FOUND AND USE_GRAPHPARTITIONING )

message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "Graph Partitioning"
                       "" )                   
	# Metis
    scai_summary_message ( "FOUND"
                           "METIS_FOUND"
                           "OPTIONAL"
                           "Metis"
                           "Version ${METIS_VERSION} at ${METIS_INCLUDE_DIR}" )

	# ParMetis
    scai_summary_message ( "FOUND"
                           "PARMETIS_FOUND"
                           "OPTIONAL"
                           "ParMetis"
                           "Version ${PARMETIS_VERSION} at ${PARMETIS_INCLUDE_DIR}" )

# EXAMPLES
message ( STATUS "" )
scai_summary_message ( "USE"
                       "BUILD_EXAMPLES"
                       "OPTIONAL"
                       "Examples"
                       "" )

# LAMA TEST
set ( REQUIRED_FOUND FALSE )
if    ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )
  set ( REQUIRED_FOUND TRUE )
endif ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )

message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "Test"
                       "" )

    # Boost Test-Framework
    scai_summary_message ( "FOUND"
                           "Boost_UNIT_TEST_FRAMEWORK_FOUND"
                           "OPTIONAL"
                           "Boost Unit Test"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
                           
    # Boost Regex
    scai_summary_message ( "FOUND"
                           "Boost_REGEX_FOUND"
                           "OPTIONAL"
                           "Boost Regex"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

# DOC
set ( REQUIRED_FOUND FALSE )
if    ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )
  set ( REQUIRED_FOUND TRUE )
endif ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )
     
message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "Documentation"
                       "" )
    # Sphinx                               
    scai_summary_message ( "FOUND"
                           "SPHINX_FOUND"
                           "OPTIONAL"
                           "Sphinx"
                           "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )

    # DOXYGEN
    scai_summary_message( "FOUND"
                          "DOXYGEN_FOUND"
                          "OPTIONAL"
                          "Doxygen"
                          "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )

message ( STATUS "" )

scai_status_message ( HEADLINE "Configuration Details:" )
message ( STATUS "" )
message ( STATUS "LAMA (ALL) Version : ${SCAI_LAMA_ALL_VERSION} ${SCAI_VERSION_NAME}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "Library Type : ${SCAI_LIBRARY_TYPE}" )
message ( STATUS "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
message ( STATUS "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" )
message ( STATUS "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
if    ( USE_CODE_COVERAGE )
	message ( STATUS "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
message ( STATUS "" )
