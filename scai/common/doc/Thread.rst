.. _Thread:

Thread
======

For threading we use an own class. In a previous version we used the thread class of boost 
that caused a lot of problems during installation as differnt libraries are
needed for differnt Boost versions. The thread class provided by the
common library uses either pThreads or the
std::thread class of C++11.

 * ThreadId used for logging
 * Thread-private variables for tracing
 * condtion, mutex for the implementations of a thread pool
 * threads that run one fucntion for asynchronous exections

The thread class provived by the common library contains only a subset of the C++11 standard.
The functionality is mainly needed for the SCAI tasking library. 
