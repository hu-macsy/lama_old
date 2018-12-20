.. _Distribution:

Distributions
=============

A distribution object is used to map data structures like vectors or matrices to the different processors
of a distributed system. It stands for a mapping of N elements (range 0 to N-1, also called global indexes)
to the processors of a corresponding communicator, so that each processor has a number of local indexes for 
which it is the owner.

The class *Distribution* itself is an abstract base class with many pure methods that must be implemented by
the derived classes, especially operations to get the local indexes on a processor and to determine the owners
for a set of global indexes.

Distributons are always created on the heap and managed by shared pointers. Therefore different distributed
data structures can share the mapping and the lifetime of the distribution ends with the lifetime of the last
object that uses it.

.. code-block:: c++

   typedef std::shared_ptr<Distribution> DistributionPtr;

   DistributionPtr dist( new BlockDistribution( N, comm ) );
   DenseVector<double> v1( dist );
   DenseVector<double> v2( dist );
   DenseMatrix<double> m( dist, dist );

.. _dmemo-distributions:

Distribution Classes
--------------------

LAMA provides the following derived distribution classes.

Block Distribution
^^^^^^^^^^^^^^^^^^

The *BlockDistribution* creates contiguous blocks of the same size (except from the last block), which are successively
assigned to the processors.

.. figure:: _images/block_distribution.*
    :width: 500px
    :align: center
  
    Block distribution of 11 elements onto 3 processors (block size is 4).
    
A BlockDistribution is created by passing the global distribution size and a communicator (optional, 
default is the current communicator). Beside the constructor a function is provided that creates
the shared pointer object.

.. code-block:: c++

   DistributionPtr dist( new BlockDistribution( N, comm ) );
   auto dist = std::make_shared<BlockDistribution>( N, comm );
   auto dist = blockDistributon( N );

Cyclic Distribution
^^^^^^^^^^^^^^^^^^^

The *CyclicDistribution* creates stripes of the given chunk size and assigns them consistently.

.. figure:: _images/cyclic2_distribution.* 
    :width: 500px
    :align: center
    :alt: CyclicDistribution

    Cyclic(2) distribution of 11 elements onto 3 processors.

You create a CyclicDistribution with the shown chunk size of '2' this way:
    
.. code-block:: c++

   DistributionPtr cyclic( new CyclicDistribution( N, 2, comm ) );

General Block Distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^

The *GenBlockDistribution* create blocks of the given sizes and assigns them successivly to the processes. The sum of
the sizes have to match the global size.

The following example creates with three parts of size 1, 3 and 2 rows/columns:

.. code-block:: c++

   HArray<IndexType> sizes( { 3, 5, 3 } );
   DistributionPtr genBlock( new GenBlockDistribution( N, sizes, comm ) );

.. figure:: _images/genblock_distribution.* 
    :width: 500px
    :align: center
    :alt: GenBlockDistribution

    General block distribution of 11 elements onto 3 processors with sizes (3, 5, 3)

Beside this constructor it is also possible to create a general block distribution by the local size
or by a weight. The following example shows how to set individual weights for the processors
by an environment variable and to use this weight for some kind of load distribution.

.. code-block:: c++

    const IndexType N = 1000;
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    float weight = 0.5;

    common::Settings::setRank( comm->getRank() );
    common::Settings::getEnvironment( weight, "SCAI_WEIGHT" );

    auto dist = genBlockDistributionByWeight( N, weight, comm );

.. code-block:: c++

    export SCAI_WEIGHT=0.3,0.5,0.7
    mpirun -np 3 ...

By these values processor 0 will have 200 values, processor 1 333 and processor 2 467 values.

GeneralDistribution
^^^^^^^^^^^^^^^^^^^

With the *GeneralDistribution* a fully free distribution can be created. 
Therefore, an array with the owner for each global index must be specified.
The array must only be available on one processor (root).

.. code-block:: c++

   HArray<PartitionId> owners( { 1, 2, 0, 1, 0, 0, 2, 2, 1, 1, 1 } );  // 11 entries
   PartitionId root = 0;
   DistributionPtr gen = generalDistributionByOwners( owners, root );
   
In this example process 0 owns index 2, 4, and 5, process 1 owns 0, 3, 8, 9, and 10 and 
process 2 owns 1, 6, 7.

An alternative constructor uses the individual sets of owned indexes on each processor.

.. code-block:: c++

    IndexType N = 11;

    HArray<IndexType> myIndexes;

    switch ( comm->getRank() ):
    {
        case 0 : myIndexes = HArray<IndexType>( { 2, 4, 5 } ); 
                 break;
        case 1 : myIndexes = HArray<IndexType>( { 0, 3, 8, 9, 10 } );
                 break;
        case 2 : myIndexes = HArray<IndexType>( { 1, 6, 7 } );
                 break;
    }

    auto gen2 = generalDistribution>( N, myIndexes, comm );

For the latter constructor the number of locally owned indexes must sum up to the global size and
each global index must appear exactly once in the local array ``myIndexes`` on a processor. It is not possible
that one element is owned by multiple processors.

.. figure:: _images/general_distribution.* 
    :width: 500px
    :align: center
    :alt: GeneralDistribution

    General distribution of 11 elements onto 3 processors.

Compared to the other distributions, general distributions have the big disadvantage that one local processor
does not know the full mapping, i.e. it cannot determine the owner of an abritrary index without further
communication. This also implies that each method requiring the computation of ownership must be
called by all processors.

Grid Distribution
^^^^^^^^^^^^^^^^^

A *GridDistribution* stands for a block distribution of an n-dimenisonal grid in multiple dimensions.

.. code-block:: c++

    const IndexType N1 = 5;
    const IndexType N2 = 4;
    Grid globalGrid( N1, N2 );
    Grid procGrid( 2, 2 );
    DistributionPtr gridDist( new GridDistribution( globalGrid, comm, procGrid ) );

Actually, this defines a mapping from the indexes 0 to N1 * N2 - 1 to four processors. The elements of the
grid are assumed to be stored in a row-major order, i.e. ( x, y+1 ) follows directly ( x, y ) and 
there are N2 elements between ( x + 1, y ) and ( x, y ).

.. figure:: _images/grid_distribution.* 
    :width: 700px
    :align: center
    :alt: GridDistribution

The number of processors in the processor grid has to match the size of the communicator, i.e. the number
of processors onto which the application is running. The procGrid argument is optional in the constructor
of a grid distribution. If it is not specified a processor grid is built from the available processors
in such a way that an optimal balancing with smallest boundaries is achieved.

Single Distribution
^^^^^^^^^^^^^^^^^^^

A *SingleDistribution* stands for a mapping of a all data to one single processor, i.e. only one
processor owns all the data.
    
.. code-block:: c++

    const PartitionId p = 2;
    DistributionPtr singleDist( new SingleDistribution( p, comm ) );

Joined Distribution
^^^^^^^^^^^^^^^^^^^

A *JoinedDistribuiton* is the concatenation of two mappings.
    
.. code-block:: c++

    auto dist1 = blockDistribution( N1 );
    auto dist2 = blockDistribution( N2 );
    auto dist = joinedDistributions( dist1, dist2 );   // mapping for N1 + N2 elements

Even if it stands on its own for a distribution, it becomes especially useful for joined data structures
where the joined data is not built explicitly and exists only implicitly.

No Distribution
^^^^^^^^^^^^^^^

Since there are cases you need to assign a *DistributionPtr* to a constructor or function, but you do not want to
distribute the data (in one direction) you have the possibility to create a *NoDistribution*. It invokes that there is
no distribution of the data and all processes have a local copy.

.. code-block:: c++

   DistributionPtr no( new NoDistribution ( numRows ) );

Regarding distributed memory programming you should keep in mind that not distributed data might either be used
in a private mode where each processor works on individual values or in a global mode, where all processors
have exactly the same values for their incarnation.

Methods for Distributions
-------------------------

In the following some important methods of a distribution are shortly described and explained.
For a detailed description of the virtual methods of a distribution we refer to the reference documentation.

Owned Indexes
^^^^^^^^^^^^^

All distributions provides a method to get an array with all global indexes that are owned by
the corresponding processor. This method does never required any communication and can be called
at any time.

.. code-block:: c++

    DistributionPtr dist = ...;
    HArray<IndexType> ownedIndexes;
    dist->getOwnedIndexes( ownedIndexes );

Computation of Ownership
^^^^^^^^^^^^^^^^^^^^^^^^

For operations on distributed data structures it might be the case that elements from other processors
are required, i.e. elements that reside on another processors. One important step for the
corresponding communications is to compute the owners of these required elements. Therefore each distribution
implements a method to query the owners for an array of global indexes.

.. code-block:: c++

    HArray<IndexType> requiredIndexes = ...;
    DistributionPtr dist = ...
    HArray<PartitionId> owners;
    dist->computeOwners( owners, requiredIndexes );

The following figure shows a typical example of such a call. Each processor calls this
method with its individual set of requried indexes to get the owners.

.. figure:: _images/compute_owners.* 
    :width: 500px
    :align: center
    :alt: ComputeOwnersBlock

    Computation of ownership with a block distribution of 40 elements onto 4 processors.

While for most distributions it is a simple operation to compute the ownership, e.g. for a
block distribution it is just an  integer divide operation, it can be rather complex for
general distribution where it also involves communication.

One possible solution is to build on each a processor an array that contains the owner for
each global index. While this is the most efficient solution it has the big disadvantage that
it might require too much memory, especially for a very large number of processors.

Another less efficient solution is to set up a block distributed array of all owners. Each
processor asks for its required indexes the corresponding processors for the owners.
This corresponds a gather operation of a distributed array.

.. figure:: _images/compute_owners_general.* 
    :width: 700px
    :align: center
    :alt: ComputeOwnersGeneral

    Computation of ownership with a general distribution of 40 elements onto 4 processors.

Comparison of Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Usually, many data structures will be distributed among the available processors, e.g. two vectors might be distributed.
For the implementation of operations on these distributed data structures, it is important to know whether two data
structures have the same distribution, as in such a case many operations can be implemented without any 
communication at all.

.. code-block:: c++

   DistributionPtr d1( new GenBlockDistribution ( n ) );
   DistributionPtr d2( new GenBlockDistribution ( n ) );

   ...

   if ( *d1 == *d2 )
   {
      // implement the operation on the local parts
      ....
   } 
   else
   {
       COMMON_THROWEXCEPTION( "Operation not available, different distributions" )
   }
