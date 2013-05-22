###
 # @file CMakeLists.txt
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
 # @brief lama/CMakeLists.txt
 # @author Jan Ecker
 # @date 28.03.2013
 # $Id$
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        BaseAccess
        CommunicationPlan
        Communicator
        CommunicatorFactory
        CommunicatorManager
        Context
        ContextAccess
        ContextFactory
        ContextManager
        DefaultHostContext
        DefaultHostContextManager
        DenseVector
        Distributed
        DistributionFactory
        DistributionManager
        HostContext
        LAMAArray
        LAMAArrayUtils
        LAMAArrayView
        LAMAInterface
        LAMAInterfaceRegistry
        NoCommunicator
        NoCommunicatorManager
        Printable
        ReadAccess
        Scalar
        StorageIO
        SyncToken
        Vector
        Walltime
        WriteAccess
    )

set ( HEADERS
        config
        HostReadAccess
        HostWriteAccess
        LAMATypes
        NonCopyable
        NoSyncToken
        tracing
        TypeTraits
    )

lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )

### add subdirectories ###
if ( CUDA_FOUND AND LAMA_USE_CUDA )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND LAMA_USE_CUDA )
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
add_subdirectory ( norm )
add_subdirectory ( openmp )
add_subdirectory ( solver )
add_subdirectory ( storage )
add_subdirectory ( task )
add_subdirectory ( tracing )

### set link directories ###
if ( WIN32 )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

link_directories ( ${LAMAPP_LINK_DIRECTORIES} )

### add library ###
add_library ( lama SHARED ${CXX_SOURCES} ${HEADERS} )

### set target link libraries ###
set( LAMA_LINK_LIBRARIES log4lama )

if ( CUDA_FOUND AND LAMA_USE_CUDA )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} lamacuda )
endif ( CUDA_FOUND AND LAMA_USE_CUDA )

if ( MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} lamampi )
endif ( MPI_FOUND AND LAMA_USE_MPI )

if ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_PBLAS_LIBRARIES} ${LAMA_BLAS_LIBRARIES})
    add_definitions( -DLAMA_MKL_SCALAPACK )
else ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_BLAS_LIBRARIES} )
endif ( SCALAPACK_FOUND AND MPI_FOUND AND LAMA_USE_MPI )

if ( NOT WIN32 )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${Boost_THREAD_LIBRARY} )
endif ( NOT WIN32 )

if ( LAMA_BLAS_FOUND )
    set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${LAMAPP_LINK_LIBRARIES} )
endif( LAMA_BLAS_FOUND )

# add additional link libraries
set ( LAMA_LINK_LIBRARIES ${LAMA_LINK_LIBRARIES} ${LAMA_ADDITIONAL_LINK_LIBRARIES} )

target_link_libraries ( lama ${LAMA_LINK_LIBRARIES} )

if ( LAMA_ADDITIONAL_LINK_FLAGS )
    set_target_properties ( lama PROPERTIES LINK_FLAGS ${LAMA_ADDITIONAL_LINK_FLAGS} )
endif ( LAMA_ADDITIONAL_LINK_FLAGS )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DLAMA_COMPILING_DLL )
endif ( WIN32 )

### install ###
install ( TARGETS lama DESTINATION lib )
