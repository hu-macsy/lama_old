###
 # @file scai/lama/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief Main specification file for CMake
 # @author Thomas Brandes, Jiri Kraus
 # @date 29.03.2012
###

### Set PROJECT NAME / supported LANGUAGES
project ( scai_lama )

### includes ###

## important Variables for project
set ( SCAI_LAMA_VERSION       2.0.0 )
include ( SetProjectVariables)
set ( SCAI_LAMA_EXTERNAL_DEPS OpenMP ZLIB PNG ) # 11
include ( Dependencies/minimalSupportedVersions )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
include ( SetBuildFlags )
include ( VersionCheck/Compiler )
include ( Compiler/SetCPPFlags )

## Use Boost if BUILD_TEST is enabled
if ( BUILD_TEST )
    include( Package/BoostTest )
endif ( BUILD_TEST )

### include all external packages ( optional )

include( Package/OpenMP )
include( Package/ZLIB )
include( Package/PNG )

message( STATUS "OpenMP: ${SCAI_OPENMP_INCLUDE_DIR} ${SCAI_OPENMP_LIBRARIES}" )
message( STATUS "ZLIB: ZLIB_FOUND=${ZLIB_FOUND}, USE_ZLIB=${USE_ZLIB}, ${SCAI_ZLIB_INCLUDE_DIR} ${SCAI_ZLIB_LIBRARIES}" )
message( STATUS "PNG: ${SCAI_PNG_INCLUDE_DIR} ${SCAI_PNG_LIBRARIES}" )

## finalize build flags
include ( ConcludeFlags )

### add definitions ###

include ( Settings/logging )
include ( Settings/tracing )

add_definitions ( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )
add_definitions ( -DSCAI_TRACING_${SCAI_TRACING} )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include ( also configured files )
set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )

include_directories ( ${SCAI_PNG_INCLUDE_DIR} )
include_directories ( ${SCAI_ZLIB_INCLUDE_DIR} )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        _DenseVector
        _SparseVector

        DenseVector
        GridVector
        GridSection
        SparseVector

        Vector

        GridVector
        GridSection

    HEADERS                  # .hpp only

        Scalar
        GridReadAccess
        GridWriteAccess
    )

### add subdirectories
add_subdirectory ( doc )

add_subdirectory ( expression )
add_subdirectory ( io )
add_subdirectory ( matrix )
add_subdirectory ( matutils )
add_subdirectory ( mepr )
add_subdirectory ( norm )
add_subdirectory ( storage )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

### set link directories ###
if ( WIN32 )
	link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${SCAI_LAMA_VERSION} )

## add link libraries

target_link_libraries( ${PROJECT_NAME} scai_common scai_logging scai_tracing scai_tasking scai_hmemo 
                                       scai_kregistry scai_blaskernel scai_utilskernel scai_sparsekernel scai_dmemo )

# external libraries

target_link_libraries ( ${PROJECT_NAME} ${SCAI_OPENMP_LIBRARIES} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_ZLIB_LIBRARIES} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_PNG_LIBRARIES} )

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION ${SCAI_LAMA_INCLUDE_DIR} )

## install library.hpp
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/lama.hpp.in" "${CONFIGURE_INCLUDE_DIR}/scai/lama.hpp" )
install ( FILES "${CONFIGURE_INCLUDE_DIR}/scai/lama.hpp" DESTINATION include/scai )

