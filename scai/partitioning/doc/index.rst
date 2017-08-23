.. _main-page_partitioning:

#################
SCAI Partitioning
#################

***********
Description
***********

The **Partitioning** is a library that provides interfaces to graph and mesh partitioning
software like Metis, Parmetis, and Zoltan.

* A partioning is a class that provides a routine to compute the row distribution for a sparse matrix
  in such a way that a good load balancing is achieved and that minimizes halo communication.

********
Contents
********

Here is a list of provided classes of the Partitioning library

======================== ================================================================================
Class                    Description
======================== ================================================================================
:ref:`Partitioning`      Abstract base class for partioning
:ref:`BlockPartitioning` Simple partitioning only for load distribution
:ref:`MetisPartitioning` Partitioning using the Metis software.
======================== ================================================================================

.. toctree::
   :hidden:

   Partitioning
   BlockPartitioning
   MetisPartitioning

*******
Example
*******

Here is a short example:

.. code-block:: c++

    #include <scai/partitoning.hpp>

    using namespace scai;
    using namespace lama;
    using namespace partitioning;

    ...

    CSRSparseMatrix<double> A( fileName );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitioningPtr partitioning( Partitioning::create( "METIS" ) );
    DistributionPtr dist( partitioning->partitionIt( comm, A, 1.0f ) );

    A.redistribute( dist, dist );

*********************
Environment Variables
*********************

None.

************
Dependencies
************

Internal dependencies:

* :ref:`SCAI Common<scaicommon:main-page_common>`
* :ref:`SCAI Dmemo<scaidmemo:main-page_dmemo>`

External dependencies: 

* :ref:`Metis`

.. toctree::
   :hidden:

   Metis
