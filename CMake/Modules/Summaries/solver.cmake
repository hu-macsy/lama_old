###
 # @file Summaries/solver.cmake
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

include ( Functions/scaiStatusMessage )
include ( Functions/scaiSummaryMessage )

include ( VersionDefinition )
include ( CompilerVersion )
include ( CheckC++11 )

message ( STATUS "" )
message ( STATUS "Summary of solver Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

# C++ Compiler
scai_status_message ( HEADLINE "Compiler:" )

scai_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )
                       
message ( STATUS "" )

if    ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )
    set( REQUIRED_FOUND TRUE )
else  ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )
	set( REQUIRED_FOUND FALSE )
endif ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "Solver"
                       "Needs compiler supporting C++11 or Boost and pThreads" )

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

# LAMA (core)
message ( STATUS "" )
scai_status_message ( HEADLINE "LIBRARIES:" )
   
#    # Boost
#    scai_summary_message ( "FOUND"
#                           "BOOST_INCLUDE_DIR"
#                           "Boost"
#                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${BOOST_INCLUDE_DIR}" )

# LAMA TEST
message ( STATUS "" )
scai_summary_message ( "USE"
                       "BUILD_TEST"
                       "LAMA TEST"
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

# LAMA CUDA
message ( STATUS "" )
scai_summary_message ( "FOUND"
                       "SCAI_COMMON_FOUND"
                       "SCAI common"
                       "" )
                       
scai_summary_message ( "FOUND"
                       "SCAI_LOGGING_FOUND"
                       "SCAI logging"
                       "" )
                       
scai_summary_message ( "FOUND"
                       "SCAI_TRACING_FOUND"
                       "SCAI tracing"
                       "" )

scai_summary_message ( "FOUND"
                       "SCAI_TASKING_FOUND"
                       "SCAI tasking"
                       "" )
                       
scai_summary_message ( "FOUND"
                       "SCAI_KREGISTRY_FOUND"
                       "SCAI kregistry"
                       "" )

scai_summary_message ( "FOUND"
                       "SCAI_BLASREGISTRY_FOUND"
                       "SCAI blasregistry"
                       "" )
                           
scai_summary_message ( "FOUND"
                       "SCAI_HMEMO_FOUND"
                       "SCAI hmemo"
                       "" )

scai_summary_message ( "FOUND"
                       "SCAI_LAMA_FOUND"
                       "SCAI lama"
                       "" )

message ( STATUS "" )
scai_status_message ( HEADLINE "DOCUMENTATION:" )
# DOC
scai_summary_message ( "USE"
                       "BUILD_DOC"
                       "DOC"
                       "" )
                                     
scai_summary_message ( "FOUND"
                       "SPHINX_FOUND"
                       "Sphinx"
                       "Version ${Sphinx_VERSION_STRING} at ${Sphinx-build_EXECUTABLE}: 'make doc' to build user documentation" )

message ( STATUS "" )

message ( STATUS "" )

scai_status_message ( HEADLINE "INFO:" )
message ( STATUS "LAMA Version : ${LAMA_VERSION} ${LAMA_VERSION_NAME}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "Library Type : ${SCAI_LIBRARY_TYPE}" )
message ( STATUS "ASSERT Level : ${SCAI_ASSERT_LEVEL}" )
message ( STATUS "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" )
message ( STATUS "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
if    ( USE_CODE_COVERAGE )
	message ( STATUS "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
message ( STATUS "" )


# Check if all required packages are found
# LAMA (core)

if    ( NOT ( SCAI_BLAS_FOUND AND BOOST_INCLUDE_DIR  ) OR ( (SCAI_BLAS_NAME MATCHES "BLAS") AND NOT LAPACK_FOUND )  )
    message( FATAL_ERROR "Configuration for LAMA (core) incomplete!")
endif ( NOT ( SCAI_BLAS_FOUND AND BOOST_INCLUDE_DIR  ) OR ( (SCAI_BLAS_NAME MATCHES "BLAS") AND NOT LAPACK_FOUND ) )

# LAMA MPI
if    ( USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "Configuration for LAMA MPI incomplete!")
endif ( USE_MPI AND NOT MPI_FOUND )

# LAMA MPI
if    ( USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "Build of LAMA MPI enabled, but configuration is incomplete!")
endif ( USE_MPI AND NOT MPI_FOUND )

# LAMA Cuda
if    ( USE_CUDA AND NOT CUDA_FOUND )
    message( FATAL_ERROR "Build of LAMA Cuda enabled, but configuration is incomplete!")
endif ( USE_CUDA AND NOT CUDA_FOUND )

# LAMA Test
if    ( LAMA_BUILD_TEST AND NOT ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND ) )
    message( FATAL_ERROR "Build of LAMA Test enabled, but configuration is incomplete!")
endif ( LAMA_BUILD_TEST AND NOT ( Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND ) )
