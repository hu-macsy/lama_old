### set include directories ### 

include_directories ( . ${CMAKE_SOURCE_DIR}/framework/src )

if ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIRS} )
else ( Boost_INCLUDE_DIR )
    message ( FATAL_ERROR "Boost_INCLUDE_DIR not found" )
endif ( Boost_INCLUDE_DIR )

### set source files ###

#WARNING: Put the *Benchmark.cpp IN FRONT of the InputSet*.cpp files!
#         Otherwise you will end up getting awkward segmentation faults.
set ( BENCHSOURCES
        Registry.cpp
        LAMAInputSet.cpp
        LAMAInputSetCreators.cpp
        LAMAInputSetComplexityVisitor.cpp
        LAMARandomInputSetCreator.cpp
        LAMAFileInputSetCreator.cpp
        ConvertMatrixStorageBenchmark.cpp
     )

# The parallel benchmarks will be built in any case but 
# should only be executed in case of MPI or PGAS interface is available

set ( BENCHSOURCES
         ${BENCHSOURCES}
         PCGBenchmark.cpp
         PGMRESBenchmark.cpp
         PJacobiBenchmark.cpp
         PSimpleAMGCGBenchmark.cpp
         PSimpleAMGComponentsBenchmark.cpp
         PSimpleAMGSetupBenchmark.cpp
         PMetaSolverBenchmark.cpp
         PMVBenchmark.cpp
         LAMAMPIBenchmark.cpp
     )

### set target link libraries ###  

set ( BENCH_LINK_LIBRARIES benchmark_framework lama log4lama )

if ( WIN32 )
    set ( LAMA_BENCH_LINK_DIRECTORIES ${LAMA_BENCH_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

link_directories ( ${LAMA_BENCH_LINK_DIRECTORIES} )

### Glob for headers to make them available in Visual Studio Projects ###

file ( GLOB_RECURSE HEADERS *.hpp )

### add definitions ###

if ( WIN32 )
    add_definitions ( -DLAMABENCH_COMPILING_DLL )
endif ( WIN32 )


### add library ###

add_library ( lama_benchmarks SHARED ${BENCHSOURCES} ${HEADERS} )

### add executable ###

add_executable ( DumpInputSet DumpInputSet.cpp ${HEADERS} )

### set target link libraries ###

target_link_libraries ( lama_benchmarks ${BENCH_LINK_LIBRARIES} )

target_link_libraries ( DumpInputSet lama_benchmarks ${BENCH_LINK_LIBRARIES} benchmark_framework )

### install ####

file ( GLOB IncludeBenchFiles *.hpp )
install ( FILES ${IncludeBenchFiles} DESTINATION bench/include )

### benchmark library should relatively be the same for CMAKE_INSTALL_PREFIX and CMAKE_BUILD_DIR

install ( TARGETS lama_benchmarks DESTINATION bench )

if ( MPI_FOUND OR ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) )
    install ( TARGETS DumpInputSet DESTINATION bin )
endif ( MPI_FOUND OR ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) )
