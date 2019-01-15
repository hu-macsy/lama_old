.. _main-page_dmemo:

##########
SCAI DMemo
##########

***********
Description
***********

DMemo stands for **Distributed Memory** and is a library that provides distribution and communication
routines for data structures using heterogeneous arrays.

* A communicator is an object that stands for a group of processes that might run on arbitrary nodes
  of the distributed memory platform. Methods are provided for data exchange between theses processes.
* A distribution defines a mapping of data (e.g. vectors, arrays) to the processes of a communicator.
* In contrary to the communication primitives provided by MPI, the
  methods provided by a communicator are more high-level routines that provide operations on arrays or
  vectors that involve communication, e.g. redistributions or halo exchange. Furthermore, they exploit
  C++ features like overloading and templates and by using the SCAI heterogeneous arrays they are aware 
  of valid instantions of the data to be communicated.
* Building subgroups of processors that execute different tasks on disjoint processor subsets is supported
  (splitting communicators).

********
Contents
********

Here is a list of provided classes of the DMemo library

=========================== ================================================================================
Class                       Description
=========================== ================================================================================
:ref:`Communicator`         Base class for communication between different partitions
:ref:`CollectiveFile`       Base class for concurrent I/O of processors
:ref:`NoCommunicator`       Default Communicator to be used on serial machines
:ref:`MPICommunicator`      MPI Communicator
:ref:`Distribution`         Mapping of an index range to a number of partitions
:ref:`CommunicationPlan`    Number of contiguous elements to exchange betweeen processors
:ref:`GlobalExchangePlan`   Communication schedule for global exchange 
:ref:`GlobalAddressingPlan` Communication schedule for gathering from and scattering into distributed data
:ref:`RedistributePlan`     Communication schedule for global redistribution of data
:ref:`HaloExchangePlan`     Communication schedule for update of halo/shadow values
=========================== ================================================================================

.. toctree::
   :hidden:

   Communicator
   CollectiveFile
   NoCommunicator
   MPICommunicator
   Distribution
   CommunicationPlan
   GlobalExchangePlan
   GlobalAddressingPlan
   RedistributePlan
   HaloExchangePlan

*******
Example
*******

Here is a short example:

.. code-block:: c++

    #include <scai/dmemo/BlockDistribution.hpp>

    using namespace scai;

    // use the default communicator
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    IndexType size = 71;
    dmemo::DistributionPtr dist ( new dmemo::BlockDistribution( size, comm ) );

Please note that objects of the classes Distribution and Communicator should
always be created as shared pointers as these objects have typically multiple owners,
e.g. a distribution can be used for multiple vectors, and a communicator for different
distributions.

*********************
Environment Variables
*********************

The default communicator is usually that communication library that has been
used for the installation. If both are supported, it can be chosen:

* ``SCAI_COMMUNICATOR`` ("MPI" for MPI parallelism or "NO" for serial execution without MPI)

When using processor arrays (2D, 3D) a default topology can be specified that
is used to split up the available processes. The total number of processors must
match the number of processors used in the default communicator.

* ``SCAI_NP`` ("2x4", "2x2x2")

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
