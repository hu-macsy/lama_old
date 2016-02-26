.. _main-page_dmemo:

##########
SCAI DMemo
##########

*************
Specification
*************

* Handles distributed memory
* Provides Communicators and Distributions
* Internal dependencies: common, logging, tracing, tasking, hmemo
* External dependencies: MPI, GPI, GraphPartitioning

********
Contents
********

.. toctree::
   :titlesonly:
   :maxdepth: 1
   
   distributions

**********
Motivation
**********

DMemo stands for **Distributed Memory** and is a library that provides distribution and communication.
routines.

*************
DMemo Classes
*************

Here is a complete list of provided classes of the DMemo library

=================     ================================================================================
Class                 Description
=================     ================================================================================
Communicator          Base class for communication between different partitions
NoCommunicator        Default communicator 
MPICommunicator       MPI Communicator
GPICommunicator       GASPI Communicator
Distribution          Mapping of an index range to a number of partitions
=================     ================================================================================

