.. _main-page_hmemo:

##########
SCAI HMemo
##########

************
Introduction 
************

HMemo stands for **Heterogeneous Memory** and is a library that provides a C++ container class
HArray (**Heterogeneous Array**) that manages multiple incarnations of the data on different devices.

Beside this container class some other classes are provided:

 * Context class to deal with different devices (Host for CPU, CUDA for GPU devices, MIC for Intel Xeon Phi, …)
 * Memory class to deal with different memory locations (device memory, CPU memory, pinned memory, …)
 * ReadAccess and Write Acess on a HArray in order to keep track at which locations valid memory data is available

**************************************
Heterogeneous Memory Library Reference
**************************************

Here is a list of all provided classes of the HMemo library

=================     ================================================================================
Class                 Description
=================     ================================================================================
:ref:`Context`        Base class for different devices
HostContext           Derived context class for working on CPU (host)
CUDAContext           Derived context class for working on CUDA devices
:ref:`Memory`         Base class for memory management at a certain location
HostMemory            Derived memory class for CPU memory management
CUDAMemory            Derived memory class for CPU memory management
CUDAHostMemory        Derived memory class for pinned host memory
HData                 One incarnation of a HArray at one memory location
HDataManager          Deals with read and write requests and initiates corresponding memory transfers
_HArray               Common base class for HArray
:ref:`HArray`         Template container class
:ref:`Access`         Read- and write access
ReadAccess            Template class for read access on HArray
WriteAccess           Template class for write access on HArray
=================     ================================================================================

.. toctree::
   :hidden:
   
   HArray
   Context
   Memory
   Access

************
Dependencies
************

Internal dependencies:

* :ref:`SCAI Common<scaicommon:main-page_common>`
* :ref:`SCAI Logging<scailogging:main-page_logging>`
* :ref:`SCAI Tracing<scaitracing:main-page_tracing>`
* :ref:`SCAI Tasking<scaitasking:main-page_tasking>`

External dependencies:

* `CUDA <http://www.nvidia.com/object/cuda_home_new.html>`_ for CUDA Context
* Compiler supporting Intel MIC Architecture for using the Xeon Phi Coprocessor
