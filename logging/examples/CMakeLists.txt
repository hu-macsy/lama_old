
find_package(OpenMP)

if (OPENMP_FOUND)

    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

endif()

add_executable( LogTest LogTest.cpp )
target_link_libraries( LogTest ${LOGGING_LIBRARY} )
set_target_properties( LogTest PROPERTIES COMPILE_FLAGS -DLAMA_LOG_LEVEL_TRACE )

add_executable( ThreadTest ThreadTest.cpp )
target_link_libraries( ThreadTest ${LOGGING_LIBRARY} )

if ( OPENMP_FOUND )
    # Adding this example only possible if OpenMP is available
    add_executable( OpenMPTest OpenMPTest.cpp )
    target_link_libraries( OpenMPTest ${LOGGING_LIBRARY} )
endif()

