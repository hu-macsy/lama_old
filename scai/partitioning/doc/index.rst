.. _main-page_partitioning:

#################
SCAI Partitioning
#################

***********
Description
***********

The **Partitioning** is a library that provides interfaces to graph and mesh partitioning
software like Metis, Parmetis, and Zoltan.

* A partioning is a class that provides routines to compute the distribution for a sparse matrix
  in such a way that a good load balancing is achieved and that minimizes either the size 
  of the halo storage or the amount of communication for the exchange of halo data, i.e. non-local
  values required by the individual processors.

********
Contents
********

Here is a list of provided classes of the Partitioning library

============================ ================================================================================
Class                        Description
============================ ================================================================================
:ref:`Partitioning`          Abstract base class for partioning
:ref:`BlockPartitioning`     Simple block partitioning only for load distribution
:ref:`CyclicPartitioning`    Simple cyclic partitioning only for load distribution
:ref:`MetisPartitioning`     Partitioning using the software package Metis
:ref:`ParMetisPartitioning`  Partitioning using the software package ParMetis.
============================ ================================================================================

.. toctree::
   :hidden:

   Partitioning
   BlockPartitioning
   CyclicPartitioning
   MetisPartitioning
   ParMetisPartitioning

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
