###
 # @file scai/dmemo/CMakeLists.txt
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
 # @brief CMake specification file for project dmemo (distributed memory)
 # @author Thomas Brandes
 # @date 16.02.2016
###

### Set PROJECT NAME / supported LANGUAGES
project ( scai_dmemo )

set ( INTERNAL_DEPS common logging tracing tasking hmemo kregistry blaskernel utilskernel sparsekernel )
set ( SCAI_DMEMO_EXTERNAL_DEPS OpenMP MPI GPI GraphPartitioning ) 

string ( REGEX REPLACE "([a-z]+)" "${SCAI_LIBRARY_PREFIX}\\1" INTERNAL_LIBS "${INTERNAL_DEPS}")

### includes ###

## important Variables for project
include ( SetProjectVariables)
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
include( Package/MPI )
include( Package/GPI )
include( Package/GraphPartitioning )

message( STATUS "OpenMP: ${SCAI_OPENMP_INLCUDE_DIR} ${SCAI_OPENMP_LIBRARIES}" )
message( STATUS "MPI: ${SCAI_MPI_INLCUDE_DIR} ${SCAI_MPI_LIBRARIES}" )
message( STATUS "GPI: ${SCAI_GPI_INLCUDE_DIR} ${SCAI_GPI_LIBRARIES}" )
message( STATUS "GraphPartitioning: ${SCAI_GRAPHPARTITIONING_INLCUDE_DIR} ${SCAI_GRAPHPARTITIONING_LIBRARIES}" )

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

include_directories ( ${SCAI_OPENMP_INCLUDE_DIR} )
include_directories ( ${SCAI_MPI_INCLUDE_DIR} )
include_directories ( ${SCAI_GPI_INCLUDE_DIR} )
include_directories ( ${SCAI_GRAPHPARTITIONING_INCLUDE_DIR} )

if    ( ( METIS_FOUND OR PARMETIS_FOUND ) AND USE_GRAPHPARTITIONING )
    if    ( SCAI_INDEX_TYPE STREQUAL "int" )
        set ( METISDISTRIBUTION MetisDistribution )
    endif ( SCAI_INDEX_TYPE STREQUAL "int" )
endif ( ( METIS_FOUND OR PARMETIS_FOUND ) AND USE_GRAPHPARTITIONING )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        Communicator
        NoCommunicator
        Distributed
        CommunicationPlan

        BlockDistribution
        CyclicDistribution
        Distribution
        GenBlockDistribution
        GeneralDistribution
        GridDistribution
        Halo
        HaloBuilder
        NoDistribution
        Redistributor
        SingleDistribution

        ${METISDISTRIBUTION}

    )

### add subdirectories

add_subdirectory ( doc )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

if ( GPI_FOUND AND USE_GPI )
    include( Package/Boost )
    message( STATUS "add BOOST_INCLUDE_DIR=${BOOST_INCLUDE_DIR}" )
    include_directories ( ${BOOST_INCLUDE_DIR} )
    add_subdirectory( gpi )
endif ( GPI_FOUND AND USE_GPI )

if ( MPI_FOUND AND USE_MPI )
    add_subdirectory ( mpi )
endif ( MPI_FOUND AND USE_MPI )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### set link directories ###
if ( WIN32 )
	link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${SCAI_VERSION} )

## add link libraries

target_link_libraries ( ${PROJECT_NAME} ${INTERNAL_LIBS} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_OPENMP_LIBRARIES} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_MPI_LIBRARIES} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_GPI_LIBRARIES} )
target_link_libraries ( ${PROJECT_NAME} ${SCAI_GRAPHPARTITIONING_LIBRARIES} )

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install "<library>.hpp"
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

## add custum target distClean
include ( CustomCommands/DistClean )

