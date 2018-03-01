.. _Partitioning:

Partitioning
============

Graph partitioning is a well studied problem in combinatorial scientific computing. 
An important application is the distribution of rows of a sparse matrix for a
matrix vector multiplication where the goals are to balance the load and to minimize 
communication.

LAMA itself does not implement own graph partitioning algorithms but provides
interfaces to well established software packages like Metis.

Each (derived) partioning class must implement two pure methods, one for partitioning
of square symmetric matrices and one for arbirtrary rectangular matrices.

.. code-block:: c++

    virtual void squarePartitioning( hmemo::HArray<PartitionId>& mapping,
                                     const lama::_Matrix& matrix,
                                     const hmemo::HArray<float>& processorWeights ) const = 0;

    virtual void rectangularPartitioning( hmemo::HArray<PartitionId>& rowMapping,
                                          hmemo::HArray<PartitionId>& colMapping,
                                          const lama::_Matrix& matrix,
                                          const hmemo::HArray<float>& processorWeights ) const = 0;

The array ``processorWeights`` contains a weight factor for each target processor. The number of
target processors does not have to match the number of processors onto which the input matrix
is distributed. So both methods might be called on a single processor to compute a mapping for an
arbitrary number of processors.

The input matrix must always be a sparse matrix, in the first method it must also be square.
It might be distributed among the available processors on which the application runs. Each processor
returns for its local part the new mapping, i.e. how the rows (and columns for the rectangular partitioning)
will be mapped to the specfified number of processors (size of processorWeights).
