
enable_language ( C )

set ( CMAKE_THREAD_PREFER_PTHREAD 1 )
set ( THREADS_PREFER_PTHREAD_FLAG 1 )

find_package( Threads ) # use ${CMAKE_THREAD_LIBS_INIT} for target_link_libraries

###  Here we use PThread library for threads
###  Note: FindThreads in CMake is available as Module, but is buggy, needs update of CheckIncludeFiles.cmake
#find_library ( PTHREADS_LIBRARY NAMES pthread pthreads )