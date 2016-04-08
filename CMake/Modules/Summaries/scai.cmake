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
include ( Settings/bashFormats )

emptyline()
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )

heading ( "Compiler:" )

if    ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND TRUE )
else  ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND FALSE )
endif ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )

heading2 ( "Configuration" "REQUIRED_FOUND" )

found_message ( "C++ Compiler" "CMAKE_CXX_COMPILER" "REQUIRED" "" )
found_message ( "with C++11 support" "CXX_SUPPORTS_C11" "REQUIRED" "" )

if    ( NOT CXX_SUPPORTS_C11 )
    emptyline()
    message ( STATUS "Either compiler supporting C++11 or Boost needed." )
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )

# LAMA (core)
heading ( "External Libraries:" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )
    set ( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )

heading2 ( "Required core" "REQUIRED_FOUND" )

    found_message ( "pThreads" "SCAI_THREAD_LIBRARIES" "REQUIRED" "Version ${SCAI_THREAD_VERSION}" )
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
    found_message ( "BLAS" "SCAI_BLAS_FOUND" "REQUIRED" "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )
    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
        message ( STATUS "                                 ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        found_message ( "Lapack" "LAPACK_FOUND" "REQUIRED" "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )

heading2 ( "Optional components" "" )
#indent_message ( "2" "Optional components" )
#emptyline()

# OpenMP usage

heading3 ( "OpenMP" "USE_OPENMP" )

    found_message ( "OpenMP" "OPENMP_VERSION" "OPTIONAL" "Version ${OPENMP_VERSION}" )
    found_message ( "compile flag" "OpenMP_CXX_FLAGS" "OPTIONAL" "${OpenMP_CXX_FLAGS}" )
    found_message ( "schedule type" "SCAI_OMP_SCHEDULE" "OPTIONAL" "set to \"${SCAI_OMP_SCHEDULE}\"" )

# LAMA CUDA
set ( REQUIRED_FOUND FALSE )
if    ( CUDA_FOUND AND USE_CUDA )
  set ( REQUIRED_FOUND TRUE )
endif ( CUDA_FOUND AND USE_CUDA )

heading3 ( "CUDA" "USE_CUDA" )

    found_message ( "CUDA" "CUDA_FOUND" "OPTIONAL" "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
    found_message ( "Compute Capability" "CUDA_HAVE_GPU" "OPTIONAL" "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
heading3 ( "MIC" "USE_MIC" )

# LAMA MPI
set ( REQUIRED_FOUND FALSE )
if    ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )
  set ( REQUIRED_FOUND TRUE )
endif ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )

heading3 ( "Graph Distributed" "REQUIRED_FOUND" )

    found_message ( "MPI" "MPI_FOUND" "OPTIONAL" "Version ${MPI_VERSION} at ${SCAI_MPI_INCLUDE_DIR}" )
    found_message ( "GPI" "GPI_FOUND" "OPTIONAL" "at ${SCAI_GPI_INCLUDE_DIR}" )

# Graph Partitioning
set ( REQUIRED_FOUND FALSE )
if    ( METIS_FOUND AND USE_GRAPHPARTITIONING )
  set ( REQUIRED_FOUND TRUE )
endif ( METIS_FOUND AND USE_GRAPHPARTITIONING )

heading3 ( "Graph Partitioning" "REQUIRED_FOUND" )

    found_message ( "Metis" "METIS_FOUND" "OPTIONAL" "Version ${METIS_VERSION} at ${METIS_INCLUDE_DIR}" )
    found_message ( "ParMetis" "PARMETIS_FOUND" "OPTIONAL" "Version ${PARMETIS_VERSION} at ${PARMETIS_INCLUDE_DIR}" )

# EXAMPLES
heading3 ( "Examples" "BUILD_EXAMPLES" )

# LAMA TEST
set ( REQUIRED_FOUND FALSE )
if    ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )
  set ( REQUIRED_FOUND TRUE )
endif ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )

heading3 ( "Test" "REQUIRED_FOUND" )

    found_message ( "Boost Unit Test" "Boost_UNIT_TEST_FRAMEWORK_FOUND" "OPTIONAL" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
    found_message ( "Boost Regex" "Boost_REGEX_FOUND" "OPTIONAL" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

# DOC
set ( REQUIRED_FOUND FALSE )
if    ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )
  set ( REQUIRED_FOUND TRUE )
endif ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )

heading3 ( "Documentation" "REQUIRED_FOUND" )

    found_message ( "Sphinx" "SPHINX_FOUND" "OPTIONAL" "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )
    found_message ( "Doxygen" "DOXYGEN_FOUND" "OPTIONAL" "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )

heading ( "Configuration Details:" )

indent_message ( "1" "LAMA (ALL) Version : ${SCAI_LAMA_ALL_VERSION} ${SCAI_VERSION_NAME}" )
indent_message ( "1" "Build Type   : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Library Type : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
indent_message ( "1" "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" )
indent_message ( "1" "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
if    ( USE_CODE_COVERAGE )
    indent_message ( "1" "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
emptyline()
