###
 # @file scai/common/CMakeLists.txt
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
 # @brief CMake configuration file for common library
 # @author Thomas Brandes
 # @date 09.06.2015
 # @since 2.0.0
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
    include ( Compiler/VersionCheck )
    include ( Compiler/SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )
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
include ( SetCacheVariableTexts )

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

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...
scai_project (

    CLASSES                  # .cpp, .hpp

        Printable
        ScalarType
        SCAITypes
        ContextType
        Settings
        Thread
        Walltime
        LibModule
        TypeTraits

    HEADERS                  # .hpp only

        Complex
        config
        Constants
        Factory
        Factory1
        Math
        OpenMP
        preprocessor
        ReductionOp
        NonCopyable
)

## configure c++11/boost dependent sources
set ( CONFIGURE_DIR "${CONFIGURE_INCLUDE_DIR}/scai/${PROJECT_SURNAME}" )

if    ( NOT CXX_SUPPORTS_C11 )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/boost_check.hpp.in" "${CONFIGURE_DIR}/boost_check.hpp" )
    install ( FILES ${CONFIGURE_DIR}/boost_check.hpp DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

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
install ( FILES ${CONFIGURE_DIR}/bind.hpp ${CONFIGURE_DIR}/function.hpp
                ${CONFIGURE_DIR}/shared_ptr.hpp ${CONFIGURE_DIR}/unique_ptr.hpp
                ${CONFIGURE_DIR}/weak_ptr.hpp DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

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
