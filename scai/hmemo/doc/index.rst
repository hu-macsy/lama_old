.. _main-page_hmemo:

##########
SCAI HMemo
##########

*************
Specification 
*************

* Handles memory
* Supports different backends (currently Host, CUDA, MIC)
* Provides the Context and HArray
* Internal dependencies: common, logging, tracing, tasking

**********
Motivation
**********

HMemo stands for **Heterogeneous Memory** and is a library that provides a C++ container class
HArray (**Heterogeneous Array**) that manages multiple incarnations of the data on different devices.

Beside this container class some other classes are provided:

 * Context class to deal with different devices (Host for CPU, CUDA for GPU devices, MIC for Intel Xeon Phi, …)
 * Memory class to deal with different memory locations (device memory, CPU memory, pinned memory, …)
 * ReadAccess and Write Acess on a HArray in order to keep track at which locations valid memory data is available

*************
HMemo Classes
*************

Here is a complete list of all provided classes of the HMemo library

=================     ================================================================================
Class                 Description
=================     ================================================================================
Context               Base class for different devices
HostContext           Derived context class for working on CPU (host)
CUDAContext           Derived context class for working on CUDA devices
Memory                Base class for memory management at a certain location
HostMemory            Derived memory class for CPU memory management
CUDAMemory            Derived memory class for CPU memory management
CUDAHostMemory        Derived memory class for pinned host memory
HData                 One incarnation of a HArray at one memory location
HDataManager          Deals with read and write requests and initiates corresponding memory transfers
_HArray               Common base class for HArray
HArray                Template container class
ReadAccess            Template class for read access on HArray
WriteAccess           Template class for write access on HArray
=================     ================================================================================

.. toctree::
   :titlesonly:
   :maxdepth: 2
   
   HArray
   Context
   Memory
   Access



