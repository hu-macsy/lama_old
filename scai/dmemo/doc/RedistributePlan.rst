.. _RedistributePlan:

RedistributePlan
================

A redistribution of data is a very common operation for distributed memory programming,
i.e. a distributed array with a source distribution becomes assigned to a distributed
array with a target distribution. 

.. code-block:: c++

    HArray<ValueType> sourceArray = ...;   // local part of distributed array via sourceDist
    HArray<ValueType> targetArray;         // becomes local part of distributed array via targetDist
 
    globalAssign( targetArray, targetDist, sourceArray, sourceDist );

Actually it is straightforward to implement this operation by using a global addressing plan as
already shown. In contrary to a GlobalAddressingPlan we split up the self communication, i.e. the
data that remains local.

.. figure:: _images/redistribute_plan.*
    :width: 700px
    :align: center
    :alt: RedistributePlan

    RedistributePlan as GlobalAddressingPlan with extracted self communication.

The indexes for local communication are stored in separate arrays. Even if this is the main
distinction to a global addressing plan, there are some other differences:

 * a redistribute plan stands always for a full permutation and so it might be used 
   vice versa
 * there is an optimized constructor for building a plan from an existing source distribution
   and just new owners, i.e. the new target distribution is also built during the construction
 * the RedistributePlan supports also redistribution of sparse and dense matrices

Construction of a RedistributePlan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code shows how to build a redistribute plan from a global addressing plan.

.. code-block:: c++

    // build global plan for scattering from source dist into target dist
    auto globalPlan = globalAddressingPlan( *targetDist, sourceDist->ownedGlobalIndexes() );

    CommunicationPlan sourceCommPlan, targetCommPlan;
    HArray<IndexType> sourcePerm, targetPerm;

    globalPlan.splitUp( sourcePerm, sourceCommPlan, targetCommPlan, targetPlan );

    RedistributePlan redistPlan( targetDist, std::move( targetPerm ), std::move( targetCommPlan ),
                                 sourceDist, std::move( sourcePerm ), std::move( sourceCommPlan ) );

As this is a very common construction, a corresponding function is provided:

.. code-block:: c++

    auto plan = redistributePlanByNewDistribution( targetDist, sourceDist );

Very often a new target distribution is given by the new owner for each owned global index 
that a processor owns in the current source distribution. 

.. code-block:: c++

    HArray<PartitionId> newOwners = ...;   // specifiy for each element of data where it goes

    auto targetDist = generalDistribution( sourceDist, newOwners );
    auto plan = redistributePlanByNewDistribution( targetDist, sourceDist );

This solution is rather inefficient as it queries the owners of the owned global indexes of
the current source distribution in the target distribution. But they are actually known
and so the following constructor function is more efficient.

.. code-block:: c++

    HArray<PartitionId> newOwners = ...;   // specifiy for each element of data where it goes
    auto plan = redistributePlanByNewOwners( newOwners, dist );

The array with the new owners is used like in the global exchange plan to set up a global
exchange plan. This global exchange plan is used to send each processor its new global indexes
that it will own for the new distribution. The new global indexes coming from the different
processors form the new target distribution and will be sorted in ascending order. The
corresponding permutation can be used to unpack the redistributed data.

.. figure:: _images/redistribution.*
    :width: 700px
    :align: center
    :alt: Redistribution

    Definition of a new distribution by new owners for the local indexes of the current one.

This solution requires a sorting of all incoming global indexes from other processors to build
the array of owned global indexes for the new target distribution. This sorting might be
optimized by some kind of mergesort, as each incoming bundle from other processors is already
sorted.

Redistribution
^^^^^^^^^^^^^^

The redistribute plan can be applied to redistribute multiple vectors and/or matrices that have 
the corresponding source distribution and have afterwards the new target distribution.

.. code-block:: c++

    plan.redistribute( targetArray, sourceArray );

A dense matrix that is row distributed can be redistributed with the following method:

.. code-block:: c++

    plan.redistributeN( targetMatrix, sourceMatrix, N );

Here sourceMatrix contains the local data regarding the source distribution and targetMatrix will contain
the local data regarding the target distribution. N is the number of columns.

Very often sparse matrix data is stored like ragged arrays, especially for the wide-spread
CSR format. The row-wise distributed data has a different size for each row. The following method
redistributes such data where it is assumed that the sizes for each row is already redistributed.
 
.. code-block:: c++

    plan.redistributeV( targetArray, targetSizes, sourceArray, sourceSizes );


