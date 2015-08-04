###
 # @file Functions.cmake
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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

## Need to be macros not functions, because modifications of the parent scope

# sets the LAMA_SOURCE_DIR (used to mark the path of the actual build target
macro    ( lama_set_source_dir )
    set ( LAMA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
    set ( CXX_SOURCES "" )
endmacro ( lama_set_source_dir )

# Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures
# the installation of the header files
macro    ( lama_classes )
    lama_sources ( ${ARGN} )
    lama_headers ( ${ARGN} )
endmacro ( lama_classes )

# Adds a list of sources to the target (the related *.cpp files)
macro    ( lama_sources )    
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    
    foreach   ( SOURCE_FILE ${ARGN} )
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.cpp" )
    endforeach ( SOURCE_FILE ${ARGN} )
endmacro ( lama_sources )

# Adds a list of headers to the target configures the installation of the header files
macro    ( lama_headers )
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    get_relative_path ( RELATIVE_PATH )
    
    # clear CXX_HEADERS
    set ( CXX_HEADERS "" )
    
    foreach    ( SOURCE_FILE ${ARGN} )
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.hpp" )
        set ( CXX_HEADERS ${CXX_HEADERS} "${SOURCE_FILE}.hpp" )
    endforeach ( SOURCE_FILE ${ARGN} )
    
    # install CXX_HEADERS
    install ( FILES ${CXX_HEADERS} DESTINATION "include/${RELATIVE_PATH}" )
endmacro ( lama_headers )

# Publishes sources and headers in the parent scope
macro    ( lama_add )
    set ( CXX_SOURCES ${CXX_SOURCES} PARENT_SCOPE )
    set ( CXX_HEADERS ${CXX_HEADERS} PARENT_SCOPE )
endmacro ( lama_add )