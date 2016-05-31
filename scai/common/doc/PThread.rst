.. _PThread:

PThread Library
---------------

|Posix Threads|, usually referred to as Pthreads, is a parallel execution model, for which
nearly all Unix systems provide an API defined by the POSIX standard.

.. |Posix Threads| raw:: html

  <a href="http://en.wikipedia.org/wiki/POSIX_Threads" target="_blank"> Posix Threads </a>

PThreads can be used in C or C++ with a ``pthread.h`` header file and a pthread library.

Using PThread in LAMA
^^^^^^^^^^^^^^^^^^^^^

The common library provides a Thread class that is either implemented by using
Posix threads or by using the C++11 threads.

So a phtread library is only required if the C++ compiler does not support C++11 features.

Note: The former release LAMA 1.0 used the Boost thread library instead of Pthreads but caused
much more problems during the installation. 

Configuration with CMake
^^^^^^^^^^^^^^^^^^^^^^^^

* LAMA provides a CMAKE module file if a project uses PThreads as external library: ``Package/Thread.cmake``
