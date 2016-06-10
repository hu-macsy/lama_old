.. _PThread:

PThread Library
---------------

|Posix Threads|, usually referred to as PThreads, is a parallel execution model, for which nearly all Unix systems provide an API defined by the POSIX standard.

.. |Posix Threads| raw:: html

  <a href="http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/pthread.h.html" target="_blank"> Posix Threads </a>

PThreads can be used in C or C++ with the ``pthread.h`` header file and the pthread library.

.. Using PThread in LAMA
.. ^^^^^^^^^^^^^^^^^^^^^

The common library provides a Thread class that is implemented by using Posix threads and used by the constitutive libraries, espacially for asynchronous tasking.

.. or by using the C++11 threads.
.. So a pthread library is only required if the C++ compiler does not support C++11 features.
.. Note: The former release LAMA 1.0 used the Boost thread library instead of Pthreads but caused much more problems during the installation. 

Configuration with CMake
^^^^^^^^^^^^^^^^^^^^^^^^

The general CMake Package ``FindThreads`` does not use any hints for directing to special paths, so you should have it in your system path.

.. * LAMA provides a CMAKE module file if a project uses PThreads as external library: ``Package/Thread.cmake``
