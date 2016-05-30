###
 # @file scai/common/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration file for common library
 # @author Thomas Brandes
 # @date 09.06.2015
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_common CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### includes ###

## important Variables for project
include ( VersionDefinition )
include ( SetProjectVariables)
include ( Dependencies/internal )
include ( Dependencies/external )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( VersionCheck/Compiler )
    include ( Compiler/SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Package/doc )
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## C++11 or Boost
include ( Compiler/CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

### include all search packages
addInternalAndExternalPackages()

### finalize build flags
include ( ConcludeFlags )

### add definitions
add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )
addInternalAndExternalIncludes()
if    ( NOT CXX_SUPPORTS_C11 )
    include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 )

set ( CONFIGURE_DIR "${CONFIGURE_INCLUDE_DIR}/scai/${PROJECT_SURNAME}" )

## begin for configuring SCAITypes.hpp

set ( EXT_HOST_TYPES_LIST ${SCAI_HOST_TYPES_LIST} ) #SCAI_HOST_TYPES_LIST - long double, ComplexLongDouble
list ( REMOVE_ITEM EXT_HOST_TYPES_LIST ${LONG_DOUBLE_VALUES} )
set ( CUDA_TYPES_LIST ${EXT_HOST_TYPES_LIST} )      #SCAI_HOST_TYPES_LIST - long double, ComplexLongDouble
set ( MIC_TYPES_LIST ${EXT_HOST_TYPES_LIST} )       #SCAI_HOST_TYPES_LIST - long double, ComplexLongDouble

list ( LENGTH SCAI_HOST_TYPES_LIST NUM_ALL_HOST_TYPES )
list ( LENGTH EXT_HOST_TYPES_LIST NUM_EXT_HOST_TYPES )
list ( LENGTH CUDA_TYPES_LIST NUM_CUDA_TYPES )
list ( LENGTH MIC_TYPES_LIST NUM_MIC_TYPES )

math ( EXPR NUM_ALL_HOST_ARRAY_TYPES "${NUM_ALL_HOST_TYPES}+1" )
math ( EXPR NUM_EXT_HOST_ARRAY_TYPES "${NUM_EXT_HOST_TYPES}+1" )
math ( EXPR NUM_CUDA_ARRAY_TYPES     "${NUM_CUDA_TYPES}+1" )
math ( EXPR NUM_MIC_ARRAY_TYPES      "${NUM_MIC_TYPES}+1" )

# make strings out of list
include ( Functions/listToString )
listToString ( ", " "${SCAI_HOST_TYPES_LIST}" ALL_HOST_TYPES )
listToString ( ", " "${EXT_HOST_TYPES_LIST}" EXT_HOST_TYPES )
listToString ( ", " "${CUDA_TYPES_LIST}" CUDA_TYPES )
listToString ( ", " "${MIC_TYPES_LIST}" MIC_TYPES )

# set "biggest" common type as ScalarRepType
# last found element in TYPE_LIST_VALUES

set ( INDEX1 -1 )
set ( INDEX2 -1 )
foreach    ( ITEM ${SCAI_HOST_TYPES_LIST} )
    list ( FIND NONE_COMPLEX_VALUES ${ITEM} FOUND_INDEX1 )
    if    ( FOUND_INDEX1 GREATER INDEX1 )
        set ( INDEX1 ${FOUND_INDEX1})
    endif ( FOUND_INDEX1 GREATER INDEX1 )

    if    ( ${USE_COMPLEX} )
        list ( FIND COMPLEX_VALUES ${ITEM} FOUND_INDEX2 )
        if    ( FOUND_INDEX2 GREATER INDEX2 )
            set ( INDEX2 ${FOUND_INDEX2})
        endif ( FOUND_INDEX2 GREATER INDEX2 )
    endif ( ${USE_COMPLEX} )

endforeach ( ITEM ${SCAI_HOST_TYPES_LIST} )

if    ( ${USE_COMPLEX} )
    if    ( INDEX1 GREATER INDEX2 )
        list ( LENGTH NONE_COMPLEX_VALUES INDEX )
        math ( EXPR INDEX "${INDEX}+${INDEX1}")
        list ( GET TYPE_LIST_VALUES ${INDEX} SCAI_SCALAR_REP_TYPE ) # take biggest combinatation as value type
    else  ( INDEX1 GREATER INDEX2 )
        list ( GET COMPLEX_VALUES ${INDEX2} SCAI_SCALAR_REP_TYPE ) # take biggest complex value type
    endif ( INDEX1 GREATER INDEX2 )
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/SCAITypes.hpp.in_withComplex" "${CONFIGURE_DIR}/SCAITypes.hpp" )
else  ( ${USE_COMPLEX} )
    list ( GET NONE_COMPLEX_VALUES ${INDEX1} SCAI_SCALAR_REP_TYPE ) # take biggest none complex value type
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/SCAITypes.hpp.in" "${CONFIGURE_DIR}/SCAITypes.hpp" )    
endif ( ${USE_COMPLEX} )


## end for configuring SCAITypes.hpp


### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

set ( COMMON_SOURCES 
        SCAITypes            # SCAITypes.hpp is handled via configuration
)

### COMMON classes with .hpp and .cpp file

set ( COMMON_CLASSES

        Printable
        ScalarType
        ContextType
        Settings
        Walltime
        LibModule
        TypeTraits
)

#  Header only file of the COMMON project

set ( COMMON_HEADERS 
        ${COMPLEX_HEADER}
        config
        Constants
        Factory
        Factory1
        Math
        OpenMP
        NonCopyable
)

if    ( ${USE_COMPLEX} )
    list ( APPEND COMMON_HEADERS Complex )
endif ( ${USE_COMPLEX} )

## configure c++11/boost dependent sources

if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/boost_check.hpp.in" "${CONFIGURE_DIR}/boost_check.hpp" )
    install ( FILES ${CONFIGURE_DIR}/boost_check.hpp DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

if    ( NOT CXX_SUPPORTS_C11 )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/bind.hpp.in_boost" "${CONFIGURE_DIR}/bind.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/function.hpp.in_boost" "${CONFIGURE_DIR}/function.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/shared_ptr.hpp.in_boost" "${CONFIGURE_DIR}/shared_ptr.hpp" )
    if    ( ${Boost_VERSION} GREATER "105700" )
        configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/unique_ptr.hpp.in_boost_157" "${CONFIGURE_DIR}/unique_ptr.hpp" )
    else  ( ${Boost_VERSION} GREATER "105700" )
        configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/unique_ptr.hpp.in_boost" "${CONFIGURE_DIR}/unique_ptr.hpp" )
    endif ( ${Boost_VERSION} GREATER "105700" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/weak_ptr.hpp.in_boost" "${CONFIGURE_DIR}/weak_ptr.hpp" )
else  ( NOT CXX_SUPPORTS_C11 )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/bind.hpp.in" "${CONFIGURE_DIR}/bind.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/function.hpp.in" "${CONFIGURE_DIR}/function.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/shared_ptr.hpp.in" "${CONFIGURE_DIR}/shared_ptr.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/unique_ptr.hpp.in" "${CONFIGURE_DIR}/unique_ptr.hpp" )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/weak_ptr.hpp.in" "${CONFIGURE_DIR}/weak_ptr.hpp" )
endif ( NOT CXX_SUPPORTS_C11 )

## configure c++11/pthread dependent headers, add specific sources

if    ( CXX_SUPPORTS_C11 )
    list ( APPEND COMMON_SOURCES Thread )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/Thread.hpp.in" "${CONFIGURE_DIR}/Thread.hpp" )
else  ( CXX_SUPPORTS_C11 )
    list ( APPEND COMMON_SOURCES Thread_pthread )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/Thread.hpp.in_pthread" "${CONFIGURE_DIR}/Thread.hpp" )
endif (  CXX_SUPPORTS_C11 )

install ( FILES ${CONFIGURE_DIR}/bind.hpp ${CONFIGURE_DIR}/function.hpp
                ${CONFIGURE_DIR}/shared_ptr.hpp ${CONFIGURE_DIR}/unique_ptr.hpp
                ${CONFIGURE_DIR}/weak_ptr.hpp ${CONFIGURE_DIR}/SCAITypes.hpp
                ${CONFIGURE_DIR}/Thread.hpp 
          DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

scai_project( 

    SOURCES ${COMMON_SOURCES}
    CLASSES ${COMMON_CLASSES}
    HEADERS ${COMMON_HEADERS}
)

## add subdirectories

add_subdirectory ( cuda )
add_subdirectory ( doc )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

add_subdirectory ( exception )
add_subdirectory ( macros )
add_subdirectory ( mepr )
add_subdirectory ( mic )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## add link libraries
addExternalLinkLibraries()

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install library.hpp
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

### add custum target distClean
include ( CustomCommands/DistClean )

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/${PROJECT_SURNAME} )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
