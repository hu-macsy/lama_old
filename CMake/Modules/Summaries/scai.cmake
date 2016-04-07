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

emptyline()
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )
emptyline()

include ( Settings/bashFormats )

#scai_status_message ( HEADLINE "Compiler:" "" )
indend_message ( "1" "${TextUnderline}Compiler:${TextUnderlineReset}" )
emptyline()

if    ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND TRUE )
else  ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )
    set( REQUIRED_FOUND FALSE )
endif ( ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR) )

scai_summary_message ( "STATIC"
                       "Configuration"
                       "2"
                       "REQUIRED_FOUND"
                       "REQUIRED"
                       "" )

# C++ Compiler
scai_summary_message ( "FOUND"
                       "C++ Compiler"
                       "3"
                       "CMAKE_CXX_COMPILER"
                       "REQUIRED"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )

scai_summary_message ( "FOUND"
					             "with C++11 support"
                       "3"
                       "CXX_SUPPORTS_C11"
                       "REQUIRED"
					             "" )

if    ( NOT CXX_SUPPORTS_C11 )
    emptyline()
    indend_message ( STATUS "Either compiler supporting C++11 or Boost needed." )

    scai_summary_message ( "FOUND"
                           "Boost"
                           "3"
                           "SCAI_BOOST_INCLUDE_DIR"
                           "REQUIRED"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )

# LAMA (core)
emptyline()
#scai_status_message ( HEADLINE "External Libraries:" )
indend_message ( "1" "${TextUnderline}External Libraries:${TextUnderlineReset}" )
emptyline()

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )
    set ( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )

scai_summary_message ( "STATIC"
                       "Required core"
                       "2"
                       "REQUIRED_FOUND"
                       "REQUIRED"
                       "" )

emptyline()

    # pthreads
    scai_summary_message ( "FOUND"
                           "pThreads"
                           "3"
                           "SCAI_THREAD_LIBRARIES"
                           "REQUIRED"
                           "VERSION ${SCAI_THREAD_VERSION}" )

    # boost
    scai_summary_message ( "FOUND"
                           "Boost"
                           "3"
                           "SCAI_BOOST_INCLUDE_DIR"
                           "REQUIRED"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

    # BLAS
    scai_summary_message ( "FOUND"
                           "BLAS"
                           "3"
                           "SCAI_BLAS_FOUND"
                           "REQUIRED"
                           "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )

    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
        message ( STATUS "                                 ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
                           
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        scai_summary_message ( "FOUND"
                               "Lapack"
                               "3"
                               "LAPACK_FOUND"
                               "REQUIRED"
                               "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )

emptyline()
indend_message ( "2" "Optional components" )
emptyline()

# OpenMP usage

scai_summary_message ( "USE"
                       "OpenMP"
                       "2"
                       "USE_OPENMP"
                       "OPTIONAL"
                       ""   )

    scai_summary_message ( "FOUND"
                           "OpenMP"
                           "3"
                           "OPENMP_VERSION"
                           "OPTIONAL"
                           "Version ${OPENMP_VERSION}" )

    scai_summary_message ( "FOUND"
                           "compile flag"
                           "3"
                           "OpenMP_CXX_FLAGS"
                           "OPTIONAL"
                           "${OpenMP_CXX_FLAGS}" )

    scai_summary_message ( "FOUND"
                           "schedule type"
                           "3"
                           "SCAI_OMP_SCHEDULE"
                           "OPTIONAL"
                           "set to \"${SCAI_OMP_SCHEDULE}\"" )

# LAMA CUDA
set ( REQUIRED_FOUND FALSE )
if    ( CUDA_FOUND AND USE_CUDA )
  set ( REQUIRED_FOUND TRUE )
endif ( CUDA_FOUND AND USE_CUDA )

emptyline()
scai_summary_message ( "USE"
                       "CUDA"
                       "2"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "" )

    # CUDA
    scai_summary_message ( "FOUND"
                           "CUDA"
                           "3"
                           "CUDA_FOUND"
                           "OPTIONAL"
                           "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
                           
    # CUDA Compute Capability
    scai_summary_message ( "FOUND"
                           "Compute Capability"
                           "3"
                           "CUDA_HAVE_GPU"
                           "OPTIONAL"
                           "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
emptyline()
scai_summary_message ( "USE"
                       "MIC"
                       "2"
                       "USE_MIC"
                       "OPTIONAL"
                       "" )

# LAMA MPI
set ( REQUIRED_FOUND FALSE )
if    ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )
  set ( REQUIRED_FOUND TRUE )
endif ( ( MPI_FOUND AND USE_MPI ) OR ( GPI_FOUND AND USE_GPI ) )

emptyline()
scai_summary_message ( "USE"
                       "Distributed"
                       "3"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "" )

    # MPI
    scai_summary_message ( "FOUND"
                           "MPI"
                           "3"
                           "MPI_FOUND"
                           "OPTIONAL"
                           "Version ${MPI_VERSION} at ${SCAI_MPI_INCLUDE_DIR}" )

    # GPI
    scai_summary_message ( "FOUND"
                           "GPI"
                           "3"
                           "GPI_FOUND"
                           "OPTIONAL"
                           "at ${SCAI_GPI_INCLUDE_DIR}" )

# Graph Partitioning
set ( REQUIRED_FOUND FALSE )
if    ( METIS_FOUND AND USE_GRAPHPARTITIONING )
  set ( REQUIRED_FOUND TRUE )
endif ( METIS_FOUND AND USE_GRAPHPARTITIONING )

emptyline()
scai_summary_message ( "USE"
                       "Graph Partitioning"
                       "2"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "" )                   
	# Metis
    scai_summary_message ( "FOUND"
                           "Metis"
                           "3"
                           "METIS_FOUND"
                           "OPTIONAL"
                           "Version ${METIS_VERSION} at ${METIS_INCLUDE_DIR}" )

	# ParMetis
    scai_summary_message ( "FOUND"
                           "ParMetis"
                           "3"
                           "PARMETIS_FOUND"
                           "OPTIONAL"
                           "Version ${PARMETIS_VERSION} at ${PARMETIS_INCLUDE_DIR}" )

# EXAMPLES
emptyline()
scai_summary_message ( "USE"
                       "Examples"
                       "2"
                       "BUILD_EXAMPLES"
                       "OPTIONAL"
                       "" )

# LAMA TEST
set ( REQUIRED_FOUND FALSE )
if    ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )
  set ( REQUIRED_FOUND TRUE )
endif ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND AND BUILD_TEST )

emptyline()
scai_summary_message ( "USE"
                       "Test"
                       "2"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "" )

    # Boost Test-Framework
    scai_summary_message ( "FOUND"
                           "Boost Unit Test"
                           "3"
                           "Boost_UNIT_TEST_FRAMEWORK_FOUND"
                           "OPTIONAL"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
                           
    # Boost Regex
    scai_summary_message ( "FOUND"
                           "Boost Regex"
                           "3"
                           "Boost_REGEX_FOUND"
                           "OPTIONAL"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

# DOC
set ( REQUIRED_FOUND FALSE )
if    ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )
  set ( REQUIRED_FOUND TRUE )
endif ( ( SPHINX_FOUND OR DOXYGEN_FOUND ) AND BUILD_DOC )
     
emptyline()
scai_summary_message ( "USE"
                       "Documentation"
                       "2"
                       "REQUIRED_FOUND"
                       "OPTIONAL"
                       "" )
    # Sphinx                               
    scai_summary_message ( "FOUND"
                           "Sphinx"
                           "3"
                           "SPHINX_FOUND"
                           "OPTIONAL"
                           "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )

    # DOXYGEN
    scai_summary_message( "FOUND"
                          "Doxygen"
                          "3"
                          "DOXYGEN_FOUND"
                          "OPTIONAL"
                          "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )

emptyline()
#scai_status_message ( HEADLINE "Configuration Details:" )
indend_message ( "1" "${TextUnderline}Configuration Details:${TextUnderlineReset}" )

emptyline()
indend_message ( "1" "LAMA (ALL) Version : ${SCAI_LAMA_ALL_VERSION} ${SCAI_VERSION_NAME}" )
indend_message ( "1" "Build Type   : ${CMAKE_BUILD_TYPE}" )
indend_message ( "1" "Library Type : ${SCAI_LIBRARY_TYPE}" )
indend_message ( "1" "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
indend_message ( "1" "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" )
indend_message ( "1" "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
if    ( USE_CODE_COVERAGE )
    indend_message ( "1" "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
emptyline()
