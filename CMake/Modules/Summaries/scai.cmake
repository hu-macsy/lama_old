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
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

scai_status_message ( HEADLINE "Compiler:" )
# C++ Compiler
scai_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )

message ( STATUS "" )

if    ( ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND TRUE )
else  ( ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR) )
  set( REQUIRED_FOUND FALSE )
endif ( ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR) )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "SCAI"
                       "Needs compiler supporting C++11 or Boost" )

scai_summary_message ( "FOUND"
					             "CXX_SUPPORTS_C11"
					             "C++11 support"
					             "" )
				
if    ( NOT CXX_SUPPORTS_C11 )
    scai_summary_message ( "FOUND"
                           "BOOST_INCLUDE_DIR"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}, add include dir ${BOOST_INCLUDE_DIR} to compile your sources" )
endif ( NOT CXX_SUPPORTS_C11 )

message ( STATUS "" )

set ( OPENMP_INFO_TEXT "OpenMP schedule type is set to \"${SCAI_OMP_SCHEDULE}\"" )

scai_summary_message ( "USE"
                       "USE_OPENMP"
                       "  OpenMP usage"
                       "${OPENMP_INFO_TEXT}"   )

# LAMA (core)
message ( STATUS "" )
scai_status_message ( HEADLINE "LIBRARIES:" )

if    ( SCAI_BLAS_FOUND )
    set( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
else  ( SCAI_BLAS_FOUND )
    set( REQUIRED_FOUND FALSE )
endif ( SCAI_BLAS_FOUND )

if    ( NOT SCAI_THREAD_LIBRARIES )
    set ( REQUIRED_FOUND FALSE )
endif ( NOT SCAI_THREAD_LIBRARIES )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "SCAI (core)"
                       "" )

    # pthreads
    scai_summary_message ( "FOUND"
                           "SCAI_THREAD_LIBRARIES"
                           "pThreads"
                           "" )

# boost
message ( STATUS "" )
scai_summary_message ( "FOUND"
                       "SCAI_BOOST_INCLUDE_DIR"
                       "Boost"
                       "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}, add include dir ${SCAI_BOOST_INCLUDE_DIR} to compile your sources" )

message ( STATUS "" )

    # BLAS
    scai_summary_message ( "FOUND"
                           "SCAI_BLAS_FOUND"
                           "BLAS"
                           "${SCAI_BLAS_NAME} with: ${SCAI_SCAI_BLAS_LIBRARIES}" )
                           
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        message ( STATUS "" )
        scai_summary_message ( "FOUND"
                               "LAPACK_FOUND"
                               "LAPACK"
                               "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )

set ( REQUIRED_FOUND FALSE )
if    ( MPI_FOUND OR GPI_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( MPI_FOUND OR GPI_FOUND )

# LAMA MPI
message ( STATUS "" )
scai_summary_message ( "USE"
                       "REQUIRED_FOUND"
                       "Distributed"
                       "" )

    # MPI
    scai_summary_message ( "FOUND"
                           "MPI_FOUND"
                           "MPI"
                           "at ${SCAI_MPI_INCLUDE_DIR}" )

    # GPI
    scai_summary_message ( "FOUND"
                           "GPI_FOUND"
                           "GPI"
                           "at ${SCAI_GPI_INCLUDE_DIR}" )

# Graph Partitioning
message ( STATUS "" )
scai_summary_message ( "USE"
                       "USE_GRAPHPARTITIONING"
                       "Graph Partitioning"
                       "" )                   
	# Metis
    scai_summary_message ( "FOUND"
                           "METIS_FOUND"
                           "Metis"
                           "at ${METIS_INCLUDE_DIR}" )

	# ParMetis
    scai_summary_message ( "FOUND"
                           "PARMETIS_FOUND"
                           "ParMetis"
                           "at ${PARMETIS_INCLUDE_DIR}" )

# LAMA CUDA
message ( STATUS "" )
scai_summary_message ( "USE"
                       "USE_CUDA"
                       "CUDA"
                       "" )

    # CUDA
    scai_summary_message ( "FOUND"
                           "CUDA_FOUND"
                           "CUDA"
                           "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIRS}" )
                           
    # CUDA Compute Capability
    scai_summary_message ( "FOUND"
                           "CUDA_HAVE_GPU"
                           "Compute Capability"
                           "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
message ( STATUS "" )
scai_summary_message ( "USE"
                       "USE_MIC"
                       "MIC"
                       "" )

# LAMA TEST
message ( STATUS "" )
scai_status_message ( HEADLINE "TESTING:" )

scai_summary_message ( "USE"
                       "BUILD_TEST"
                       "TEST"
                       "" )

    # Boost Test-Framework
    scai_summary_message ( "FOUND"
                           "Boost_UNIT_TEST_FRAMEWORK_FOUND"
                           "Boost Unit Test"
                           "" )
                           
    # Boost Regex
    scai_summary_message ( "FOUND"
                           "Boost_REGEX_FOUND"
                           "Boost Regex"
                           "" )

message ( STATUS "" )
scai_status_message ( HEADLINE "DOCUMENTATION:" )
         
# DOC
message ( STATUS "" )
scai_summary_message ( "USE"
                       "BUILD_DOC"
                       "DOC"
                       "" )
                                     
scai_summary_message ( "FOUND"
                       "SPHINX_FOUND"
                       "Sphinx"
                       "Version ${Sphinx_VERSION_STRING} at ${Sphinx-build_EXECUTABLE}: 'make doc' to build user documentation" )

# DOXYGEN
#scai_summary_message ( "FOUND"
#                       "DOXYGEN_FOUND"
#                       "DOXYGEN "
#                       "Version ${DOXYGEN_VERSION} at ${DOXYGEN_EXECUTABLE}: 'make doxygendoc' to build system documentation" )

message ( STATUS "" )

scai_status_message ( HEADLINE "INFO:" )
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
