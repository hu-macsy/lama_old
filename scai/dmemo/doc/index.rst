.. _main-page_dmemo:

##########
SCAI DMemo
##########

***********
Description
***********

DMemo stands for **Distributed Memory** and is a library that provides distribution and communication
routines for data structures using heterogeneous arrays.

* A distribution defines a mapping of data (e.g. vectors, arrays) to the processors of the distributed-memory
  platform.
* In contrary to the communication primitives provided by MPI, the
  communication routines provided here are more high-level routines that provide operations on arrays or
  vectors that involve communication, e.g. redistributions or halo exchange. Furthermore, they exploit
  C++ features like overloading and templates and by using the SCAI heterogeneous arrays they are aware 
  of valid instantions of the data to be communicated.
* Building subgroups of processors that execute different tasks on disjoint processor subsets is supported.

********
Contents
********

Here is a list of provided classes of the DMemo library

======================== ================================================================================
Class                    Description
======================== ================================================================================
:ref:`Communicator`      Base class for communication between different partitions
:ref:`NoCommunicator`    Default Communicator to be used on serial machines
:ref:`MPICommunicator`   MPI Communicator
:ref:`Distribution`      Mapping of an index range to a number of partitions
:ref:`CommunicationPlan` Communication schedule for exchanging non-local values
======================== ================================================================================

.. toctree::
   :hidden:

   Distribution
   Communicator
   CommunicationPlan
   NoCommunicator
   MPICommunicator

*******
Example
*******

Here is a short example:

.. code-block:: c++

    #include <scai/dmemo/Distribution.hpp>

    using namespace scai::dmemo;

    // use the default communicator
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType size = 71;
    float weight = 1.0;
    DistributionPtr dist ( Distribution::getDistributionPtr( "CYCLIC", comm, size, weight ) );

*********************
Environment Variables
*********************

The default communicator is usually that communication library that has been
used for the installation. If both are supported, it can be chosen:

* ``SCAI_COMMUNICATOR`` ("MPI" or "NO" for no distributed communication)

If a CUDA-Aware MPI installation is available, the following environment
variable should be set:

* ``SCAI_MPI_CUDA`` (bool value, e.g. 0, 1, default is 0)

By this way, communication on heterogeneous arrays can communicate
valid data on a GPU directly without explicit copy on the host.

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

* :ref:`MPI`
* :ref:`Metis`

.. toctree::
   :hidden:

   MPI
   Metis

************
Related Work
************

* |BOOST_MPI|
* |CUDA_MPI|

.. |CUDA_MPI| raw:: html

  <a href="https://devblogs.nvidia.com/parallelforall/introduction-cuda-aware-mpi" target="_blank">CUDA Aware MPI</a>

.. |BOOST_MPI| raw:: html

  <a href="http://www.boost.org/libs/mpi" target="_blank">Boost.MPI</a>
