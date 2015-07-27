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
 # @brief lama/CMakeLists.txt
 # @author Jan Ecker
 # @date 28.03.2013
 # @since 1.0.0
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### include directories ###

include_directories( ${common_INCLUDE_DIR}  )
include_directories( ${logging_INCLUDE_DIR} )
include_directories( ${tracing_INCLUDE_DIR} )
include_directories( ${tasking_INCLUDE_DIR} )
include_directories( ${memory_INCLUDE_DIR} )
include_directories( ${cudamem_INCLUDE_DIR} )

if ( METIS_FOUND )
    include_directories( ${METIS_INCLUDE_DIR} ) 
endif ( METIS_FOUND )

if ( METIS_FOUND AND PARMETIS_FOUND )
    include_directories( ${PARMETIS_INCLUDE_DIR} ) 
endif ( METIS_FOUND AND PARMETIS_FOUND )

if ( LAMA_BLAS_NAME STREQUAL "MKL" )
    include_directories ( ${MKL_INCLUDE_DIRS} )
elseif ( NOT LAMA_BLAS_NAME STREQUAL "MKL" AND NOT INTERNALBLAS_FOUND )
    include_directories( cblas/include )
    add_subdirectory( cblas )
endif ( LAMA_BLAS_NAME STREQUAL "MKL" )

### Add classes and header files ###
set ( CLASSES
        CommunicationPlan
        Communicator
        Complex
        DenseVector
        Distributed
        LAMAArrayUtils
#        LAMAArrayView
        LAMAInterface
        LAMAInterfaceRegistry
        NoCommunicator
        Scalar
        Settings
        StorageIO
        Vector
    )

set ( HEADERS
        LAMATypes
        TypeTraits
        BLASInterface
        UtilsInterface
        CRTPCommunicator
    )

lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )


### add subdirectories ###

if ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )

if ( LAMA_USE_MIC )
    add_subdirectory ( mic )
endif ( LAMA_USE_MIC )

add_subdirectory ( distribution )
add_subdirectory ( exception )
add_subdirectory ( expression )
add_subdirectory ( io )
add_subdirectory ( macros )
add_subdirectory ( matrix )
add_subdirectory ( matutils )
if ( MPI_FOUND AND LAMA_USE_MPI )
    add_subdirectory ( mpi )
endif ( MPI_FOUND AND LAMA_USE_MPI )
#if ( ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) AND LAMA_USE_PGAS )
#    add_subdirectory( pgas )
#endif ( ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) AND LAMA_USE_PGAS )
if ( GPI2_FOUND AND IBVERBS_FOUND AND LAMA_USE_GPI2 )
    add_subdirectory( gpi )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND LAMA_USE_GPI2 )
add_subdirectory ( norm )
add_subdirectory ( openmp )
add_subdirectory ( solver )
add_subdirectory ( storage )

### set link directories ###
if ( WIN32 )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

link_directories ( ${LAMAPP_LINK_DIRECTORIES} )

### add library ###
add_library ( ama ${LAMA_LIBRARY_TYPE} ${CXX_SOURCES} ${HEADERS} )


### set target link libraries ###

set ( LAMA_USER_LINK_LIBRARIES "-L\${LAMA_ROOT}/lib -lama" )

set ( LAMA_LINK_LIBRARIES ${common_LIBRARY} )
set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lcommon" )

set ( LAMA_LINK_LIBRARIES ${logging_LIBRARY} )
set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -llogging" )

set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${tasking_LIBRARY} )
set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -ltasking" )

set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${memory_LIBRARY} )
set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lmemory" )

if ( LAMA_TRACE_LEVEL STREQUAL "ON" )
    set( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${tracing_LIBRARY} )
    set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -ltracing" )
endif ( )

# amacuda
if ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} amacuda )
    set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lamacuda" )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS}  ${CUDA_cusparse_LIBRARY}
                                                                                ${CUDA_cublas_LIBRARY} 
                                                                                ${CUDA_CUDART_LIBRARY}
                                                                                ${CUDA_CUDA_LIBRARY} )          
endif ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )

# amamic
if ( LAMA_USE_MIC )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} amamic )
    set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lamamic" )
endif ( LAMA_USE_MIC )

# amampi
if ( MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} amampi )
    set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lamampi" )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${MPI_LIBRARIES} )
endif ( MPI_FOUND AND LAMA_USE_MPI )

# amagpi
if ( GPI2_FOUND AND IBVERBS_FOUND AND LAMA_USE_GPI2 )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} amagpi )
    set ( LAMA_USER_LINK_LIBRARIES "${LAMA_USER_LINK_LIBRARIES} -lamagpi" )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${GPI2_LIBRARIES} ${IBVERBS_LIBRARIES} )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND LAMA_USE_GPI2 )

set ( LAMA_USER_LINK_FLAGS "-Wl,--whole-archive ${LAMA_USER_LINK_LIBRARIES} -Wl,--no-whole-archive" )

# blas lapack 
if ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_PBLAS_LIBRARIES} ${LAMA_BLAS_LIBRARIES})
    add_definitions( -DLAMA_MKL_SCALAPACK )
else ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_BLAS_LIBRARIES} )
endif ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )

# Metis / ParMetis
if ( METIS_FOUND )
    set( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${METIS_LIBRARY} )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${METIS_LIBRARY} )
endif ( METIS_FOUND ) 

if ( METIS_FOUND AND PARMETIS_FOUND )
    set( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${PARMETIS_LIBRARY} )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${PARMETIS_LIBRARY} )
endif ( METIS_FOUND AND PARMETIS_FOUND )

if ( LAMA_BLAS_FOUND )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${LAMAPP_LINK_LIBRARIES} )
    set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${LAMAPP_LINK_LIBRARIES} )
endif( LAMA_BLAS_FOUND )

# thread library
if ( NOT WIN32 )
    # set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${LAMA_THREAD_LIBRARY} )
    # set ( ADDITIONAL_EXTERNAL_LINKER_FLAGS ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} ${LAMA_THREAD_LIBRARY} )
endif ( NOT WIN32 )

# generate additional linker flags for the examples
foreach ( lib ${ADDITIONAL_EXTERNAL_LINKER_FLAGS} )
    get_filename_component ( libname ${lib} NAME_WE )
    get_filename_component ( libpath ${lib} PATH )
    get_filename_component ( libext ${lib} EXT )
    
    if ( NOT ${libpath} STREQUAL "" )        
        list_contains(contains "-L${libpath}" "${LAMA_USER_LINK_FLAGS}")
        if ( NOT contains )
          set ( LAMA_USER_LINK_FLAGS "${LAMA_USER_LINK_FLAGS} -L${libpath}" )
        endif ( NOT contains )
    endif ( NOT ${libpath} STREQUAL "" )
    
    if ( "${libext}" STREQUAL "" )
        set ( LAMA_USER_LINK_FLAGS "${LAMA_USER_LINK_FLAGS} -l${libname}" )
    else ( "${libext}" STREQUAL "" )
        string ( LENGTH ${libname} namelength )
        math ( EXPR namelength ${namelength}-3 )
        string ( SUBSTRING ${libname} 3 ${namelength} libname )
        set ( LAMA_USER_LINK_FLAGS "${LAMA_USER_LINK_FLAGS} -l${libname}" )
    endif ( "${libext}" STREQUAL "" )
endforeach ( lib ${LAMAPP_LINK_LIBRARIES} )

# override the linker flags if lama is build SHARED
if ( NOT ${LAMA_LIBRARY_TYPE} STREQUAL STATIC )
    set ( LAMA_USER_LINK_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_LINKER_FLAGS} -ldl ${LAMA_USER_LINK_LIBRARIES}" )
endif ( NOT ${LAMA_LIBRARY_TYPE} STREQUAL STATIC )

# generate file make.inc + install it
configure_file ( "${CMAKE_SOURCE_DIR}/examples/make.inc.in" "${PROJECT_BINARY_DIR}/examples/make.inc" )
install ( FILES "${PROJECT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

if ( LAMA_ADDITIONAL_LINK_FLAGS )
    set_target_properties ( ama PROPERTIES LINK_FLAGS ${LAMA_ADDITIONAL_LINK_FLAGS} )
endif ( LAMA_ADDITIONAL_LINK_FLAGS )


# add additional link libraries
set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${LAMA_ADDITIONAL_LINK_LIBRARIES} )

target_link_libraries ( ama ${LAMA_LINK_LIBRARIES} )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### install ###
install ( TARGETS ama DESTINATION lib )
