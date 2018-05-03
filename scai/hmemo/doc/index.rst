.. _main-page_hmemo:

##########
SCAI HMemo
##########

***********
Description 
***********

HMemo stands for **Heterogeneous Memory** and is a library that provides a C++ container class
HArray (**Heterogeneous Array**) that manages multiple incarnations of the data on different devices.

Beside this container class some other classes are provided:

 * Context class to deal with different devices (Host for CPU, CUDA for GPU devices, …)
 * Memory class to deal with different memory locations (device memory, CPU memory, pinned memory, …)
 * ReadAccess and Write Acess on a HArray in order to keep track at which locations valid memory data is available

**************************************
Heterogeneous Memory Library Reference
**************************************

Here is a list of all provided classes of the HMemo library

=================     ================================================================================
Class                 Description
=================     ================================================================================
:ref:`Context`        Polymorphic class for different devices
:ref:`Memory`         Polymorphic class for memory management at a certain location
:ref:`HArray`         Template container class
:ref:`Access`         Read- and write access
=================     ================================================================================

.. toctree::
   :hidden:
   
   HArray
   Context
   Memory
   Access

*******
Example
*******

The following example defines three heterogeneous arrays A, B, C. B and C are allocated
and defined on the host CPU, and will then be added into an array A on the GPU.

.. code-block:: c++

    const IndexType N = 10;

    HArray<double> a, b, c;

    ContextPtr host = Context::getHostPtr();
    {
        WriteOnlyAccess<double> wB( b, host, N );  // allocate on host
        WriteOnlyAccess<double> wC( c, host, N );  // allocate on host

        for ( IndexType i = 0; i < N; ++i )
        {
            wB[i] = 1.0;
            wC[i] = 2.0;
        }
    }

    ContextPtr gpu = Context::getContextPtr( Context::CUDA );
    {
        ReadAccess<double> rB( b, gpu );   // implicit allocate, transfer to GPU
        ReadAccess<double> rC( c, gpu );   // implicit allocate, transfer to GPU
        WriteOnlyAccess<double> wA( a, gpu, N );  // allocate on GPU

        double* ptrA = a.get();  // pointer to GPU data
        double* ptrB = b.get();  // pointer to GPU data
        double* ptrC = c.get();  // pointer to GPU data

        add_kernel_gpu( ptrA, ptrB, ptrC, N );  // run code on the GPU
    }

The library kregistry provides support to write this code in such a way that
the context where the addition is done can be chosen at runtime.

*****
Usage
*****

Environment variables for default context and default device.

* ``SCAI_DEVICE`` 
* ``SCAI_CONTEXT`` 

************
Dependencies
************

Internal dependencies:

* :ref:`SCAI Common<scaicommon:main-page_common>`
* :ref:`SCAI Logging<scailogging:main-page_logging>`
* :ref:`SCAI Tracing<scaitracing:main-page_tracing>`
* :ref:`SCAI Tasking<scaitasking:main-page_tasking>`

External dependencies:

* :ref:`CUDA<scaicommon:CUDA>`

************
Related Work
************

* |Hemi| is very close but restricted to Host and CUDA

.. |Hemi| raw:: html

  <a href="http://harrism.github.io/hemi/" target="_blank">Hemi</a>

