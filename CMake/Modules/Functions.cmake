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

# defined functions:
#     setAndCheckCache: Function for setting LAMA_USE_{PACKAGE_NAME} variables depending on {PACKAGE_NAME}_FOUND.
#     get_relative_path: returns the relative path to the actual directory to the CMAKE_SOURCE_DIR
#     lama_get_relative_path: returns the relative path to the actual directory to the LAMA_SOURCE_DIR (Path of the actual target)
#     lama_status_message: prints colored text messages
#     checkValue: checks whether the given value is in the value list ( pass list as "${LIST}" (doublequotes !!!) )
#     lama_generate_blanks: Simple internal helper function that generates a blank string that fits the size of an given STRING to LENGTH

# defined makros:
#     lama_set_source_dir: sets the LAMA_SOURCE_DIR (used to mark the path of the actual build target
#     lama_classes: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures the installation of the header files
#     lama_sources: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures the installation of the header files
#     lama_headers: Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures # the installation of the header files
#     lama_add: Publishes sources and headers in the parent scope
#     lama_summary_message: generates messages for lama summary page

# Function for setting LAMA_USE_{PACKAGE_NAME} variables depending on {PACKAGE_NAME}_FOUND.
# Also sets cache Variables
function ( setAndCheckCache PACKAGE_NAME )
    # Create variable names with LAMA_USE_XXX and FOUND_XXX
    set ( CACHE_VARIABLE_NAME LAMA_USE_${PACKAGE_NAME} )
    set ( FOUND_VARIABLE_NAME ${PACKAGE_NAME}_FOUND )

    # Check if cache variable is already set
    if ( DEFINED ${CACHE_VARIABLE_NAME} )
    
    # if cache variable is NOT set
    else ( DEFINED ${CACHE_VARIABLE_NAME} )
        # Check if package was found
        if ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE TRUE )
        else ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE FALSE )
        endif ( ${FOUND_VARIABLE_NAME} )
        
        # if optional parameter is set, use this one as package name for message
        if( DEFINED ARGV1 )
            set ( PACKAGE_NAME ${ARGV1} )
        endif ( DEFINED ARGV1 )
        
        # Set cache variable
        set ( ${CACHE_VARIABLE_NAME} ${USE_PACKAGE} CACHE BOOL "Enable / Disable use of ${PACKAGE_NAME}" )
    endif ( DEFINED ${CACHE_VARIABLE_NAME} )
endfunction ( setAndCheckCache )

# returns the relative path to the actual directory to the CMAKE_SOURCE_DIR
function ( get_relative_path RELATIVE_PATH )
    # get relative path
    string ( LENGTH "${CMAKE_SOURCE_DIR}" CMAKE_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
    
    if ( ${CMAKE_SOURCE_DIR_LENGTH} LESS ${CMAKE_CURRENT_SOURCE_DIR_LENGTH} )
        math ( EXPR CMAKE_SOURCE_DIR_LENGTH ${CMAKE_SOURCE_DIR_LENGTH}+1 )
        set ( PATH_SUFFIX / )
    endif ()
    
    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${CMAKE_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR_LENGTH} ${PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} ${PATH}${PATH_SUFFIX} PARENT_SCOPE )
endfunction ( get_relative_path )

# returns the relative path to the actual directory to the LAMA_SOURCE_DIR (Path of the actual target)
function ( lama_get_relative_path RELATIVE_PATH )
    # get relative path
    string ( LENGTH "${LAMA_SOURCE_DIR}" LAMA_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
   
    if ( ${LAMA_SOURCE_DIR_LENGTH} LESS ${CMAKE_CURRENT_SOURCE_DIR_LENGTH} )
        math ( EXPR LAMA_SOURCE_DIR_LENGTH ${LAMA_SOURCE_DIR_LENGTH}+1 )
        set ( PATH_SUFFIX / )
    endif ()

    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${LAMA_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${LAMA_SOURCE_DIR_LENGTH} ${PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} ${PATH}${PATH_SUFFIX} PARENT_SCOPE )
endfunction ( lama_get_relative_path )

# inspired by soci colormsg function
# prints colored text messages
function ( lama_status_message )
    string ( ASCII 27 _escape )
    # ANSI Display Atributes
    set ( ERROR "1\;31" )
    set ( WARNING "33" )
    set ( INFO "2\;32" )
    set ( HEADLINE "4" )
    
    set ( coloron FALSE )
    set ( str "" )
    foreach ( arg ${ARGV} )
        if ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
            set(str "${str}${_escape}[${${arg}}m")
            set(coloron TRUE)
        else ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
            set ( str "${str}${arg}" )
            if ( coloron )
                set ( str "${str}${_escape}[0m" )
                set ( coloron FALSE )
            endif()
            set(str "${str} ")
        endif ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
    endforeach ()
    
    message (STATUS ${str} )
endfunction( lama_status_message )


function ( checkValue SINGLEVALUE VALUELIST )
    set ( BOOLVALUE FALSE )
    foreach ( ITEM ${VALUELIST} )
        if ( ${SINGLEVALUE} MATCHES ${ITEM} )
            set ( BOOLVALUE TRUE )
        endif ( ${SINGLEVALUE} MATCHES ${ITEM} )
    endforeach( ITEM ${VALUELIST} )
    if ( NOT BOOLVALUE )
        message ( FATAL_ERROR "Selected Value ${SINGLEVALUE} is no valid choice out of ${VALUELIST}" )
    endif ( NOT BOOLVALUE )
endfunction ( checkValue SINGLEVALUE VALUELIST )

# Simple internal helper function that generates a blank string that fits the size of an given STRING to LENGTH
function ( lama_generate_blanks OUTPUT STRING LENGTH )
    string ( LENGTH "${STRING}" LAMA_STRING_LENGTH )
    # -1 for correct looping from 0 to LENGTH
    math ( EXPR LAMA_MESSAGE_BLANK_LENGTH ${LENGTH}-${LAMA_STRING_LENGTH} )
    
    set ( MESSAGE_BLANKS "")
    foreach ( LAMA_I RANGE ${LAMA_MESSAGE_BLANK_LENGTH} )
        set ( MESSAGE_BLANKS "${MESSAGE_BLANKS} " )
    endforeach( LAMA_I )
    
    set ( ${OUTPUT} ${MESSAGE_BLANKS} PARENT_SCOPE )
endfunction ( lama_generate_blanks )

## Need to be macros not functions, because modifications of the parent scope

# sets the LAMA_SOURCE_DIR (used to mark the path of the actual build target
macro ( lama_set_source_dir )
    set ( LAMA_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
    set ( CXX_SOURCES "" )
endmacro ( lama_set_source_dir )

# Adds a list of classes to the target (the related *.cpp and *.hpp files) and configures
# the installation of the header files
macro ( lama_classes )
    lama_sources ( ${ARGN} )
    lama_headers ( ${ARGN} )
endmacro ( lama_classes )

# Adds a list of sources to the target (the related *.cpp files)
macro ( lama_sources )    
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    
    foreach(SOURCE_FILE ${ARGN})
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.cpp" )
    endforeach()
endmacro ( lama_sources )

# Adds a list of headers to the target configures the installation of the header files
macro ( lama_headers )
    lama_get_relative_path ( LAMA_RELATIVE_PATH )
    get_relative_path ( RELATIVE_PATH )
    
    # clear CXX_HEADERS
    set ( CXX_HEADERS "" )
    
    foreach(SOURCE_FILE ${ARGN})
        set ( CXX_SOURCES ${CXX_SOURCES} "${LAMA_RELATIVE_PATH}${SOURCE_FILE}.hpp" )
        set ( CXX_HEADERS ${CXX_HEADERS} "${SOURCE_FILE}.hpp" )
    endforeach()
    
    # install CXX_HEADERS
    install ( FILES ${CXX_HEADERS} DESTINATION "include/${RELATIVE_PATH}" )
endmacro ( lama_headers )

# Publishes sources and headers in the parent scope
macro ( lama_add )
    set ( CXX_SOURCES ${CXX_SOURCES} PARENT_SCOPE )
    set ( CXX_HEADERS ${CXX_HEADERS} PARENT_SCOPE )
endmacro ( lama_add )

# generates messages for lama summary page
macro ( lama_summary_message MESSAGE_TYPE EXPRESSION PACKAGE_NAME ADDITIONAL_INFO )
    set ( LAMA_SUMMARY_PACKAGE_NAME_LENGTH 18 )
    if ( ${MESSAGE_TYPE} STREQUAL "FOUND" )
        set ( TYPE_TRUE "FOUND" )
        set ( TYPE_FALSE "NOT FOUND" )
        set ( TYPE_INTENT "    " )
        lama_generate_blanks ( LAMA_PACKAGE_NAME_BLANKS ${PACKAGE_NAME} ${LAMA_SUMMARY_PACKAGE_NAME_LENGTH} )
    endif ( ${MESSAGE_TYPE} STREQUAL "FOUND" )
    if ( ${MESSAGE_TYPE} STREQUAL "STATIC" )
        set ( TYPE_TRUE "REQUIRED" )
        set ( TYPE_FALSE "REQUIRED" )
        set ( TYPE_INTENT " " )
        set ( LAMA_PACKAGE_NAME_BLANKS "" )
    endif ( ${MESSAGE_TYPE} STREQUAL "STATIC" )
    if ( ${MESSAGE_TYPE} STREQUAL "USE" )
        set ( TYPE_TRUE "ENABLED" )
        set ( TYPE_FALSE "DISABLED" )
        set ( TYPE_INTENT "  " )
        set ( LAMA_PACKAGE_NAME_BLANKS "" )
    endif ( ${MESSAGE_TYPE} STREQUAL "USE" )
        
    if ( ${EXPRESSION} )
        lama_status_message ( ${TYPE_INTENT} ${PACKAGE_NAME} ${LAMA_PACKAGE_NAME_BLANKS} INFO ${TYPE_TRUE} ${ADDITIONAL_INFO})
    else ( ${EXPRESSION} )
        lama_status_message ( ${TYPE_INTENT} ${PACKAGE_NAME} ${LAMA_PACKAGE_NAME_BLANKS} ERROR ${TYPE_FALSE} )  
    endif ( ${EXPRESSION} )
endmacro ( lama_summary_message )

