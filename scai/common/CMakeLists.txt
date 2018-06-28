###
 # @file scai/common/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration file for SCAI module common
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Include required scai macros

include ( scai_macro/scai_module )
include ( scai_macro/scai_project )
include ( scai_macro/scai_subdirectories )
include ( scai_macro/scai_library )
include ( scai_function/relative_install )

## Define new SCAI module project with its internal/external dependencies

scai_module(

    MODULE_NAME   common
    EXTERNAL_DEPS OpenMP CUDA dl
)

## Common now configures some include files

set ( CONFIGURE_DIR "${CMAKE_BINARY_DIR}/include/scai/${MODULE_NAME}" )

###################################################################
## Define numeric types for which classes will be instantiated   ##
###################################################################

include ( scai_function/checkValue )
include ( scai_function/listToString )

# Possible choices for numeric types ( more than one can be selected )

set ( REAL_VALUE_TYPES    "float" "double" "long double" )
set ( COMPLEX_VALUE_TYPES "scai::ComplexFloat" "scai::ComplexDouble" "scai::ComplexLongDouble" )
set ( NUMERIC_VALUE_TYPES ${REAL_VALUE_TYPES} ${COMPLEX_VALUE_TYPES} )
set ( LONG_DOUBLE_VALUES  "long double" "scai::ComplexLongDouble" )

## define USE_COMPLEX ( might be configured )

## define SCAI_HOST_TYPES_LIST ( might be configured )

set ( SCAI_HOST_TYPES_LIST ${NUMERIC_VALUE_TYPES} CACHE STRING "Instantiation Types" )
checkValues( "${SCAI_HOST_TYPES_LIST}" "${NUMERIC_VALUE_TYPES}" )

## external host types = host types without 'long' types

set ( EXT_HOST_TYPES_LIST ${SCAI_HOST_TYPES_LIST} ) #SCAI_HOST_TYPES_LIST - long double, ComplexLongDouble
list ( REMOVE_ITEM EXT_HOST_TYPES_LIST ${LONG_DOUBLE_VALUES} )

## CUDA does not support 'long' types

set ( CUDA_TYPES_LIST ${EXT_HOST_TYPES_LIST} )

## USE_COMPLEX : set it true if at least one complex type is instantiated

set ( USE_COMPLEX FALSE )

foreach ( ITEM ${SCAI_HOST_TYPES_LIST} )
    list ( FIND COMPLEX_VALUE_TYPES ${ITEM} BOOLVALUE )
    if ( ${BOOLVALUE} GREATER -1 )
        set ( USE_COMPLEX TRUE )
    endif ()
endforeach ()

message( STATUS "USE_COMPLEX=${USE_COMPLEX}" )

# make strings out of list

listToString ( ", " "${SCAI_HOST_TYPES_LIST}" CMAKE_NUMERIC_TYPES_HOST )
listToString ( ", " "${EXT_HOST_TYPES_LIST}" CMAKE_NUMERIC_TYPES_EXT_HOST )
listToString ( ", " "${CUDA_TYPES_LIST}" CMAKE_NUMERIC_TYPES_CUDA )

# select real (non-complex) types of host types 

set ( REAL_TYPES_HOST ${SCAI_HOST_TYPES_LIST} )
foreach ( ITEM ${SCAI_HOST_TYPES_LIST} )
    list ( FIND REAL_VALUE_TYPES ${ITEM} _index )
    if ( _index LESS 0 )
        list( REMOVE_ITEM REAL_TYPES_HOST ${ITEM} )
    endif()
endforeach()

listToString ( ", " "${REAL_TYPES_HOST}" CMAKE_REAL_TYPES_HOST )

# find types that might be used for FFT

set( FFT_TYPES_HOST ${SCAI_HOST_TYPES_LIST} )
foreach ( ITEM ${SCAI_HOST_TYPES_LIST} )
    list ( FIND REAL_VALUE_TYPES ${ITEM} _index )
    if ( _index GREATER -1 )
       list( GET COMPLEX_VALUE_TYPES ${_index} _ComplexType )
    else ()
        set( _ComplexType ${ITEM} )
    endif()
    list ( FIND SCAI_HOST_TYPES_LIST ${_ComplexType} _index )
    # message ( STATUS "check type ${ITEM}, complex = ${_ComplexType}, found = ${_index}" )
    if ( _index LESS 0 )
       list ( REMOVE_ITEM FFT_TYPES_HOST ${ITEM} )
    endif()
endforeach()

listToString ( ", " "${FFT_TYPES_HOST}" CMAKE_FFT_TYPES_HOST )

###################################################################
## Define index type that is used for addressing 32 or 64 bit    ##
###################################################################

# Possible and default choices for IndexType ( only one can be selected )

scai_build_variable ( NAME      SCAI_INDEX_TYPE
                      CHOICES   "int" "long" "unsigned int" "unsigned long" 
                      DEFAULT   "int"
                      DOCSTRING "IndexType" )

# add IndexType for the array types, be careful about empty string
concatString ( "scai::IndexType" ", " "${CMAKE_NUMERIC_TYPES_HOST}"  CMAKE_ARRAY_TYPES_HOST )
concatString ( "scai::IndexType" ", " "${CMAKE_NUMERIC_TYPES_EXT_HOST}"  CMAKE_ARRAY_TYPES_EXT_HOST )
concatString ( "scai::IndexType" ", " "${CMAKE_NUMERIC_TYPES_CUDA}"  CMAKE_ARRAY_TYPES_CUDA )

# CMAKE_SCALAR_REP_TYPE : determine "biggest" common type as ScalarRepType

set ( INDEX1 -1 )
set ( INDEX2 -1 )

foreach ( ITEM ${SCAI_HOST_TYPES_LIST} )

    list ( FIND REAL_VALUE_TYPES ${ITEM} FOUND_INDEX1 )
    if    ( FOUND_INDEX1 GREATER INDEX1 )
        set ( INDEX1 ${FOUND_INDEX1})
    endif ( FOUND_INDEX1 GREATER INDEX1 )

    if    ( ${USE_COMPLEX} )
        list ( FIND COMPLEX_VALUE_TYPES ${ITEM} FOUND_INDEX2 )
        if    ( FOUND_INDEX2 GREATER INDEX2 )
            set ( INDEX2 ${FOUND_INDEX2})
        endif ( FOUND_INDEX2 GREATER INDEX2 )
    endif ( ${USE_COMPLEX} )

endforeach ()

if ( ${USE_COMPLEX} )
    if    ( INDEX1 GREATER INDEX2 )
        list ( LENGTH REAL_VALUE_TYPES INDEX )
        math ( EXPR INDEX "${INDEX}+${INDEX1}")
        list ( GET NUMERIC_VALUE_TYPES ${INDEX} CMAKE_SCALAR_REP_TYPE ) # take biggest combinatation as value type
    else  ( INDEX1 GREATER INDEX2 )
        list ( GET COMPLEX_VALUE_TYPES ${INDEX2} CMAKE_SCALAR_REP_TYPE ) # take biggest complex value type
    endif ( INDEX1 GREATER INDEX2 )
else  ()
    list ( GET REAL_VALUE_TYPES ${INDEX1} CMAKE_SCALAR_REP_TYPE ) # take biggest none complex value type
endif ()

# SCAI_USE_COMPLEX ( 0 or 1 ) is required 

if ( ${USE_COMPLEX} )
    set ( SCAI_USE_COMPLEX 1 CACHE BOOL "use complex" )
else ()
    set ( SCAI_USE_COMPLEX 0 CACHE BOOL "use complex" )
endif ()

###################################
## configuring SCAITypes.hpp     ##
################################### 

# SCAI_USE_COMPLEX is used here for conditional compilation
# SCAI_yyy_TYPES_LIST will be replaced with supported types

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/SCAITypes.hpp.in" "${CONFIGURE_DIR}/SCAITypes.hpp" )    

##############################################
## Define source file for module 'common'   ##
##############################################

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

set ( COMMON_SOURCES 
        SCAITypes            # SCAITypes.hpp is handled via configuration
)

### COMMON classes with .hpp and .cpp file

set ( COMMON_CLASSES

        AccessKind
        Printable
        ScalarType
        ContextType
        Grid
        Settings
        Walltime
        LibModule
        thread
        TypeTraits
)

#  Header only file of the COMMON project

set ( COMMON_HEADERS 
        config
        BinaryOp
        CompareOp
        Constants
        Factory
        Factory1
        Math
        MatrixOp
        OpenMP
        NonCopyable
        Stencil
        UnaryOp
        Utils
        safer_memcpy
)

if ( ${SCAI_USE_COMPLEX} )
    list ( APPEND COMMON_HEADERS Complex )
endif ( ${SCAI_USE_COMPLEX} )

## Boost used -> generate and add a boost check version file 

if ( USE_BOOST_TEST )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/boost_check.hpp.in" "${CONFIGURE_DIR}/boost_check.hpp" )
    install ( FILES ${CONFIGURE_DIR}/boost_check.hpp DESTINATION "include/scai/${MODULE_NAME}" )
endif ()

install ( FILES ${CONFIGURE_DIR}/SCAITypes.hpp
          DESTINATION "include/scai/${MODULE_NAME}" )

scai_project( 

    SOURCES ${COMMON_SOURCES}
    CLASSES ${COMMON_CLASSES}
    HEADERS ${COMMON_HEADERS}
)

## add subdirectories, add cuda even if CUDA is disabled for include file

scai_subdirectories( cuda EXAMPLES exception macros mepr TEST )

if ( CUDA_FOUND AND USE_CUDA )
    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
    set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif()

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES} )

## using std::Thread in C++ compiler with C++11 support requires linking pthread library

target_link_libraries( ${MODULE_LIBRARY} pthread )

### install ###

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )

