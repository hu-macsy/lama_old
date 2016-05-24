.. _main-page_tasking:

############
SCAI Tasking
############

***********
Description 
***********

* Handles asynchronous executions
* Supported for asynchronous computations, memory transfer and communication
* Provides the SyncToken as central concept

Tasking is a library that provides support for:

* asynchronous computations on different context devices
* asynchronous memory transfer, either in one memory class but also between different memory classes
* asynchronous communication 

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   SyncToken
   ThreadPool
   Task
   TaskSyncToken
   CUDAStreamPool
   CUDAStreamSyncToken

*******
Example
*******

.. code-block:: c++

    #include <scai/tasking/Task.hpp>

    using namespace scai;

    void f() { ... }
    void g() { ... }

    {
        tasking::Task t1( g );
        f();
        t1.wait();
    }

.. image:: _images/Tasking.png
    :align: center
    :width: 1000px

************
Dependencies
************

Internal dependencies:

* common
* logging
* tracing

************
Related Work
************

Thread pools are not provided in the C++ 11 standard, but also by Boost.
