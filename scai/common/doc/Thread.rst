.. _Thread:

Thread
======

For threading we use an own class. 
The thread class provived by the common library contains only a subset of the C++11 standard.
The functionality is mainly needed for the SCAI tasking library. 

 * ThreadId used for logging
 * Thread-private variables for tracing
 * condtion, mutex for the implementations of a thread pool
 * threads that run one function for asynchronous exections

Note: In previous versions of LAMA we used either Boost Threads or pthreads. Both of them
caused troubles on certain machines during installation. With the C++11 standard we now can
rely completely on the capabilities of the C++ compiler.
