.. _main-page_tasking:

############
SCAI Tasking
############

*************
Specification 
*************

* Handles asynchronous executions
* Supported for asynchronous computations, memory transfer and communication
* Provides the SyncToken as central concept
* Internal dependencies: common, logging, tracing

**********
Motivation
**********

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
   TaskSyncToken
   ThreadPool

