### set include directories ### 

include_directories ( . )


### set source files ###

set ( FRAMEWORKSOURCES 
        BaseInputSetRegistry.cpp
        Benchmark.cpp 
        BenchmarkPrinter.cpp
        BenchmarkRegistry.cpp 
        BFError.cpp
        BFException.cpp
        Config.cpp
        InputSet.cpp
        readDir.cpp
        frame_stdlib.cpp
	    string_helper.cpp
     )


### add executable ###

add_executable ( RunBenchmark RunBenchmark.cpp ${FRAMEWORKSOURCES} )

add_executable( RegisteredInputSets RegisteredInputSets.cpp ${FRAMEWORKSOURCES} )

add_executable( RegisteredBenchmarks RegisteredBenchmarks.cpp ${FRAMEWORKSOURCES} )


### add library ###

add_library ( benchmark_framework SHARED ${FRAMEWORKSOURCES} )


### set target link libraries ###  

set ( LIBRARIES_TO_LINK log4lama )

if ( NOT WIN32 )
        set( LIBRARIES_TO_LINK ${LIBRARIES_TO_LINK} ${CMAKE_DL_LIBS} )
endif ( NOT WIN32 )

target_link_libraries ( RunBenchmark ${LIBRARIES_TO_LINK} )

target_link_libraries ( RegisteredInputSets ${LIBRARIES_TO_LINK} )

target_link_libraries ( RegisteredBenchmarks ${LIBRARIES_TO_LINK} )

target_link_libraries( benchmark_framework log4lama )


### add definitions ###

if ( WIN32 )
	add_definitions ( -DLAMABENCHFRAME_COMPILING_DLL )
endif ( WIN32 )


### install ###

file ( GLOB IncludeBenchFrameworkFiles *.hpp )
install ( FILES ${IncludeBenchFrameworkFiles} DESTINATION bench/include/framework/src )
install ( TARGETS RunBenchmark DESTINATION bin )
install ( TARGETS RegisteredInputSets DESTINATION bin )
install ( TARGETS RegisteredBenchmarks DESTINATION bin )
install ( TARGETS benchmark_framework DESTINATION lib )
