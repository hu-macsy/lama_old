###
 # @file scai/lama/CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Thomas Brandes, Jiri Kraus
 # @date 29.03.2012
 # @since 1.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_lama CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### includes ###

## important Variables for project
include ( VersionDefinition )
include ( SetProjectVariables)
include ( InternalDependencies )
include ( ExternalDependencies )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( CompilerVersion )
    include ( SCAIAssert )
    include ( SetCPPFlags )
    include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Functions/CheckValueAtCompleteBuild )
    CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## C++11 or Boost
include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

### include all search packages
addInternalAndExternalPackages()

## finalize build flags
include ( ConcludeFlags )

### add definitions ###
add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )
add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )

## MKL FOUND
if     ( SCAI_BLAS_NAME STREQUAL "MKL" )
    add_definitions( -DSCAI_MKL_BLAS )
elseif ( INTERNALBLAS_FOUND )
    add_definitions( -DSCAI_INTERNAL_BLAS )
else   ( )
    add_definitions( -DADD_ )
endif  ( )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
include_directories ( ../.. )
addInternalAndExternalIncludes()

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        DenseVector
        HArrayUtils
        StorageIO
        Vector

    HEADERS                  # .hpp only

        UtilKernelTrait
        CSRKernelTrait
        DenseKernelTrait
        DIAKernelTrait
        ELLKernelTrait
        COOKernelTrait
        JDSKernelTrait

        LAMAKernel
        LArray
        LAMAArrayRef

        Scalar
    )

### add subdirectories
if ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
    set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

add_subdirectory ( doc )
add_subdirectory( examples )
add_subdirectory ( expression )

if ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
    add_subdirectory( gpi )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )

add_subdirectory ( io )
add_subdirectory ( matrix )
add_subdirectory ( matutils )

if ( USE_MIC )
    add_subdirectory ( mic )
endif ( USE_MIC )

add_subdirectory ( norm )
add_subdirectory ( openmp )
add_subdirectory ( storage )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### set link directories ###
if ( WIN32 )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
	link_directories ( ${LAMAPP_LINK_DIRECTORIES} )
endif ( WIN32 )

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## add link libraries
addInternalAndExternalLinkLibraries()

# external libraries
if    ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
    target_link_libraries ( ${PROJECT_NAME} ${GPI2_LIBRARIES} ${IBVERBS_LIBRARIES} )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )

if ( USE_MIC )
    set ( MKLROOT $ENV{MKLROOT} )
    set_target_properties( ${PROJECT_NAME} PROPERTIES LINK_FLAGS
        "-offload-option,mic,compiler,\"-L${MKLROOT}/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group\"" )
endif ( USE_MIC )

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install library.hpp
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/scai/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/scai/lama.hpp" DESTINATION include/scai )

# generate file make.inc + install it
#configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/examples/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" )
#install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

# copy examples in build directory + install it
#file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share )

## add custum target distClean
include ( DistClean )

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/${PROJECT_SURNAME} )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
