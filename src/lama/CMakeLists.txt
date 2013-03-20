### set include directories ### 

include_directories ( . )

if ( MKL_FOUND )
        include_directories ( ${MKL_INCLUDE_DIRS} )
endif ( MKL_FOUND )
                    

### set source files ###

set ( CXXSOURCES
        BaseAccess.cpp
        CommunicationPlan.cpp
        CommunicatorFactory.cpp
        CommunicatorManager.cpp
        Communicator.cpp
        DistributionFactory.cpp
        DistributionManager.cpp
        Context.cpp
        ContextAccess.cpp
        ContextFactory.cpp
        ContextManager.cpp
        DefaultHostContext.cpp
        DefaultHostContextManager.cpp
        DenseVector.cpp
        Distributed.cpp
        HostContext.cpp
        LAMAArray.cpp
        LAMAArrayView.cpp
        LAMAArrayUtils.cpp
        LAMAInterface.cpp
        LAMAInterfaceRegistry.cpp
        NoCommunicator.cpp
        NoCommunicatorManager.cpp
        Printable.cpp
        ReadAccess.cpp
        Scalar.cpp
        StorageIO.cpp
        SyncToken.cpp
        Vector.cpp
        WriteAccess.cpp
        # for blas library
        ${BLAS_LEVEL_CPU} ${BLAS_LEVEL} ${LAPACK}
    )

### add subdirectories ###

if ( CUDA_FOUND )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND )

# Make sure that compilation works also without MPI

add_subdirectory ( distribution )

add_subdirectory ( exception )

add_subdirectory ( matrix )

add_subdirectory ( matutils )

add_subdirectory ( norm )

add_subdirectory ( task )

# Build LAMA-MPI interface if MPI is available

if ( MPI_FOUND )
    add_subdirectory ( mpi )
endif ( MPI_FOUND )

# Build LAMA-PGAS interface if OpenShMem or GPI is available

if ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) 
    add_subdirectory ( pgas )
endif ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) 

# Make sure that compilation works also without OpenCL

if ( OPENCL_FOUND )
    add_subdirectory ( OpenCL )
endif ( OPENCL_FOUND )

add_subdirectory ( openmp )

add_subdirectory ( solver )

add_subdirectory ( storage )

add_subdirectory ( io )

# tracing: region timing and/or vampir trace calls

add_subdirectory ( tracing )


### Glob for headers to make them available in Visual Studio Projects ###
file ( GLOB_RECURSE HEADERS *.hpp )


### set link directories ###

if ( DIRECTX_FOUND )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${DIRECTX_LIBRARY_DIR} )
endif ( DIRECTX_FOUND )

if ( WIN32 )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

link_directories ( ${LAMAPP_LINK_DIRECTORIES} )


### add library ###

add_library ( lama SHARED ${CXXSOURCES} ${HEADERS} )


### set target link libraries ###  

set( LAMAPP_LINK_LIBRARIES log4lama )

if ( CUDA_FOUND )
    target_link_libraries( lama lamacuda )
endif ( CUDA_FOUND )

if ( MPI_FOUND )
    target_link_libraries ( lama lamampi )
endif ( MPI_FOUND )

if ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND )
     target_link_libraries ( lama lamapgas )
endif ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) 

if ( SCALAPACK_FOUND AND MPI_FOUND )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_PBLAS_LIBRARIES} )
    add_definitions( -DLAMA_MKL_SCALAPACK )
else ( SCALAPACK_FOUND AND MPI_FOUND )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${LAMA_BLAS_LIBRARIES} )
endif ( SCALAPACK_FOUND AND MPI_FOUND )

if ( NOT WIN32 )
    set ( LAMAPP_LINK_LIBRARIES ${LAMAPP_LINK_LIBRARIES} ${Boost_THREAD_LIBRARY} )
endif ( NOT WIN32 )

if ( LAMA_BLAS_FOUND )
    message ( STATUS "LAMA_BLAS_LIBRARIES ${LAMAPP_LINK_LIBRARIES}" )
    target_link_libraries ( lama ${LAMAPP_LINK_LIBRARIES} )
endif( LAMA_BLAS_FOUND )
           
if ( DIRECTX_FOUND )
    target_link_libraries ( lama ${DIRECTX_LIBRARIES} )
endif ( DIRECTX_FOUND )

### add definitions ###

if ( WIN32 )
    add_definitions ( -DLAMA_COMPILING_DLL )
endif ( WIN32 )

if ( MKL_FOUND )
        add_definitions ( -DLAMA_HAVE_MKL )
        # default: LAMA_FORTRAN_BLAS_STYLE_LOWERCASE
elseif ( ACML_FOUND )
        add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UPCASE )
else ()
    add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
endif ()

if ( GPI_FOUND )
    add_definitions ( -DLAMA_BUILD_GPI )
else ( GPI_FOUND )
    if ( GASNET_FOUND )
        if ( OSHMEM_FOUND )
            add_definitions ( -DLAMA_BUILD_OSHMEM )
        endif ( OSHMEM_FOUND )
    endif ( GASNET_FOUND )
endif ( GPI_FOUND )

### Installation: lib + include files ###

install ( TARGETS lama DESTINATION lib )


file ( GLOB includefiles *.hpp )
install ( FILES ${includefiles} DESTINATION include/lama )

file ( GLOB IncludeDistributeFiles mpi/*.hpp )
install ( FILES ${IncludeDistributeFiles} DESTINATION include/lama/mpi/ )

file ( GLOB IncludeCUDAFiles cuda/*.hpp )
install ( FILES ${IncludeCUDAFiles} DESTINATION include/lama/cuda )

file ( GLOB IncludePgasFiles pgas/*.hpp )
install ( FILES ${IncludePgasFiles} DESTINATION include/lama/pgas )

file ( GLOB IncludePgasFunctorFiles pgas/functor/*.hpp )
install( FILES ${IncludePgasFunctorFiles} DESTINATION include/lama/pgas/functor )

file ( GLOB IncludeDistributeFiles distribution/*.hpp )
install ( FILES ${IncludeDistributeFiles} DESTINATION include/lama/distribution )

file ( GLOB IncludeLoadBalanceFiles distribution/loadBalancing*.hpp )
install ( FILES ${IncludeLoadBalanceFiles} DESTINATION include/lama/distribution/loadBalancing )

file ( GLOB IncludeExceptionsFiles exception/*.hpp )
install ( FILES ${IncludeExceptionsFiles} DESTINATION include/lama/exception )

file ( GLOB IncludeExpressionsFiles expression/*.hpp )
install ( FILES ${IncludeExpressionsFiles} DESTINATION include/lama/expression )

file ( GLOB IncludeIOFiles io/*.hpp )
install ( FILES ${IncludeIOFiles} DESTINATION include/lama/io )

file ( GLOB IncludeMacroFiles macros/*.hpp )
install ( FILES ${IncludeMacroFiles} DESTINATION include/lama/macros )

file ( GLOB IncludeMatrixFiles matrix/*.hpp )
install ( FILES ${IncludeMatrixFiles} DESTINATION include/lama/matrix )

file ( GLOB IncludeStorageFiles storage/*.hpp )
install ( FILES ${IncludeStorageFiles} DESTINATION include/lama/storage )

file ( GLOB IncludeMatutilsFiles matutils/*.hpp )
install ( FILES ${IncludeMatutilsFiles} DESTINATION include/lama/matutils )

file ( GLOB IncludeNormFiles norm/*.hpp )
install ( FILES ${IncludeNormFiles} DESTINATION include/lama/norm )

file ( GLOB IncludeSolverFiles solver/*.hpp )
install ( FILES ${IncludeSolverFiles} DESTINATION include/lama/solver )

file ( GLOB IncludeSolverCriteriaFiles solver/criteria/*.hpp )
install ( FILES ${IncludeSolverCriteriaFiles} DESTINATION include/lama/solver/criteria )

file ( GLOB IncludeSolverLoggerFiles solver/logger/*.hpp )
install ( FILES ${IncludeSolverLoggerFiles} DESTINATION include/lama/solver/logger )

file ( GLOB IncludeSolverCreatorFiles solver/creator/*.hpp )
install ( FILES ${IncludeSolverCreatorFiles} DESTINATION include/lama/solver/creator )

file ( GLOB IncludeTaskFiles task/*.hpp )
install ( FILES ${IncludeTaskFiles} DESTINATION include/lama/task )

file ( GLOB IncludeTracingFiles tracing/*.hpp )
install ( FILES ${IncludeTracingFiles} DESTINATION include/lama/tracing )
