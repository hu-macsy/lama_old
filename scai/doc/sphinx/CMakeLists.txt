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
 # @brief Main specification file for sphinx doc
 # @author Lauretta Schubert
 # @date 11.02.2016
 # @since 2.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( sphinx_dox )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../CMake/Modules" )

include ( Package/Sphinx )
include ( VersionDefinition )
include ( SetProjectVariables )

set ( SCAI_COMMON_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-common-${SCAI_COMMON_VERSION}/ ) # 1
set ( SCAI_LOGGING_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-logging-${SCAI_LOGGING_VERSION}/ ) # 2
set ( SCAI_TRACING_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-tracing-${SCAI_TRACING_VERSION}/ ) # 3
set ( SCAI_TASKING_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-tasking-${SCAI_TASKING_VERSION}/ )# 4
set ( SCAI_HMEMO_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-hmemo-${SCAI_HMEMO_VERSION}/ ) # 5
set ( SCAI_KREGISTRY_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-kregistry-${SCAI_KREGISTRY_VERSION}/ ) # 6
set ( SCAI_BLASKERNEL_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-blaskernel-${SCAI_BLASKERNEL_VERSION}/ ) # 7
set ( SCAI_DMEMO_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-dmemo-${SCAI_DMEMO_VERSION}/ ) # 8
set ( SCAI_LAMA_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-lama-${SCAI_LAMA_VERSION}/ ) # 9
set ( SCAI_SOLVER_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/scai-solver-${SCAI_SOLVER_VERSION}/ ) # 10
	    
# Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
# must be in same directory as Sphinx source files
configure_file ( "${CMAKE_SOURCE_DIR}/conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/conf.py" )

file ( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html )

# Add custom command for building the userdoc
add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    COMMAND ${Sphinx-build_EXECUTABLE} -b ${SCAI_DOC_TYPE} -c ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
	    
# Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
add_custom_target (
    doc ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
)
	    
# Add documentation to the install target
install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
    DESTINATION ${DOC_ROOT_DIR}
)
