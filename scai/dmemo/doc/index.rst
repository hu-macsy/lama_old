.. _main-page_dmemo:

##########
SCAI DMemo
##########

* Handles distributed memory
* Provides Communicators and Distributions

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

*************
Relationships
*************

The communicator is part of a distribution as it specifies the target processors onto which
the data is distributed.

*******
Example
*******

Here is a short example:

.. code-block:: c++

    #include <scai/dmemo/Distribution.hpp>

    // use the default communicator
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType size = 71;
    float weight = 1.0;
    DistributionPtr dist ( Distribution::getDistribution( "CYCLIC", comm, size, weight ) );

*********************
Environment Variables
*********************

The default communicator is usually that communication library that has been
used for the installation. If both are supported, it can be chosen:

* ``SCAI_COMM`` (either MPI or GPI)

************
Dependencies
************

Internal dependencies:

* :ref:`SCAI Common<scaicommon:main-page_common>`
* :ref:`SCAI Logging<scailogging:main-page_logging>`
* :ref:`SCAI Tracing<scaitracing:main-page_tracing>`
* :ref:`SCAI Tasking<scaitasking:main-page_tasking>`
* :ref:`SCAI Hmemo<scaihmemo:main-page_hmemo>`

External dependencies: 

* `MPI <https://www.mpi-forum.org/docs/docs.html>`_
* `GPI <http://www.gpi-site.com/gpi2>`_
* `Metis Graph Partitioning Software <http://glaros.dtc.umn.edu/gkhome/views/metis>`_ 

************
Related Work
************

* Boost MPI 
