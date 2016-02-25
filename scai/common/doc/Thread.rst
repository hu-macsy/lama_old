Thread
======

For threading we use an own class. In a previous version we used the thread class of boost 
that caused a lot of problems during installation as differnt libraries are
needed for differnt Boost version. Our class either uses pThreads or the
std::thread class of C++11.

 * ThreadId used for logging
 * Thread-private variables
 * condtion, mutex 
 * threads that run one fucntion

