###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief CMake configuration file for scai_utilskernel tests
 # @author Eric Schricker
 # @date 18.02.2016
###

set ( CXX_SOURCES
      UtilsKernelTest
      HArrayUtilsTest
      UtilsTest
      LArrayTest
) 

include_directories ( ${Boost_INCLUDE_DIR} )

add_executable ( utilsKernelTest ${CXX_SOURCES} )

target_link_libraries ( utilsKernelTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )

#if    ( CUDA_FOUND AND USE_CUDA )
#    add_subdirectory ( cuda )
#endif ( CUDA_FOUND AND USE_CUDA )
#
#if    ( USE_MIC )
#    add_subdirectory ( mic )
#endif ( USE_MIC )

### Copy the xmltest.sh file in build directory ### 
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverage.sh" @ONLY)
endif ( USE_CODE_COVERAGE )
