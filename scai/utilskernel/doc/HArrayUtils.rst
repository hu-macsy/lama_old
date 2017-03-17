Class HArrayUtils
=================

The class HArrayUtils provides a lot of arithmetic and other operations for
heterogeneous array. All methods of this class are static methods that take
heterogenous arrays as arguments.

Arithmetic Operations
---------------------

For arithmetic operations the ValueType must be known at compile time.

.. code-block:: c++

    HArray<double> x1;
    HArray<double> x2;
    HArray<double> result;

    // elementwise addition of elements
    HArrayUtils::binaryOp( result, x1, x2, binaryOp::ADD );

    // result[i] = 1.0 / result[i]
    HArrayUtils::compute( result, 1.0, binaryOp::DIVIDE, result );

    // result[i] *= 2.5
    HArrayUtils::compute( result, result, binaryOp::MULT, 2.5 );

Filling
-------

Setting all array elements with one value.

.. code-block:: c++

    HArray<double> x;

Heterogeneous arrays can be filled with random numbers. By an additional
argument, a fill rate can be specified to generate sparse arrays.

.. code-block:: c++

    _HArray& x = ...
    float fillRate = 0.1;  
    HArrayUtils::setRandom( x, 100, fillRate )

Reductions
----------

.. code-block:: c++

    HArray<double> x;

    double sum = scan( x );
    double sum = reduce( x, binaryOp::MIN );   // get min val

    double maxDiff = HArrayUtils::reduce2( x, y, binaryOp::SUB, binaryOp::ABS_MAX ); 

Conversions
-----------

Conversion operations can take arrays of different value types.

.. code-block:: c++

    HArray<double> x;
    HArray<float> y;
    HArray<IndexType> indexes;

    // x = y with conversion of each element
    HArrayUtils::assign( x, y );

    // x[indexes[i]] = y[i] 
    HArrayUtils::scatter( x, indexes, y, binaryOp::COPY );

    // x[i] = y[indexes[i]}
    HArrayUtils::scatter( x, indexes, y, binaryOp::COPY );

Sorting
-------

The sort utility can be used to sort values in ascending or descending order.

.. code-block:: c++

    HArray<double> x( ... )
    ...
    HArrayUtils::sort( x, true );   // sorts ascending
    HArrayUtils::sort( x, false );  // sorts descending

If it is important to know the original positions of the sorted values, an additional
argument for the permutation can be specified.

.. code-block:: c++

    HArray<double> x( ... )
    HArray<double> x0( x );   // save old values
    HArray<IndexType> perm;
    ...
    HArrayUtils::sort( x, perm, true );   // sorts ascending

    // x1 = x0[ perm ], x1 has the same sorted values as x
    HArray<double> x1;
    HArrayUtils::gather( x1, x0, perm, binary::COPY )

The permutation array can be inverted.

.. code-block:: c++

    // x1[invPerm] = x0, x1 has the same sorted values as x

    HArray<IndexType> invPerm;
    HArrayUtils::setInversePerm( invPerm, perm );
    HArrayUtils::scatter( x1, x0, invPerm, binary::COPY )

Furthermore, a predicate ``isSorted`` is provided to check if the values of an array are sorted.

.. code-block:: c++

    SCAI_ASSERT_ERROR( HArrayUtils::isSorted( x, true ), "Sorting (ascending) was not correct" )

Beside this general sort routine, a more efficient bucket sort routine is provided.

.. code-block:: c++

    const IndexType numBuckets;
    HArray<IndexType> array = ...

    // array contains only values in the range 0, ..., numBuckets-1
    // e.g. numBuckets = 6, array = [ 0, 5, 2, 1, 3, 4, 4, 3, 1, 2, 5, 0 ]
  
    HArray<IndexType> offsets;   // will be the offsets
    HArray<IndexType> perm;      // is the permutation 

    HArrayUtils::bucketSort( offsets, perm, array, numBuckets, loc );

    // offsets = [0,      2,    4,    6,    8,    10,    12 ]
    // perm    = [ 0, 11, 3, 8, 2, 9, 4, 7, 5, 6, 2, 10 ]

The output array perm contains the indexes how the values have to be sorted.
The output array offsets will have the size numBuckets + 1, and for an arbirtrary bucket i
the array perm[offsets[i]], ..., perm[offsets[i+1]-1] contains the indexes of the array values
belonging to the bucket i.

Sparse Arrays
-------------

The following operations can be used to build a sparse array from a dense array
and vice versa.

.. code-block:: c++

    HArray<double> denseArray;
    HArray<IndexType> sparseIndexes;
    HArray<double> sparseArray;

    HArrayUtils::buildSparseIndexes( sparseIndexes, denseArray );

    // by gathering the non-zero values an array with the non-zero values is built

    HArrayUtils::setGather( sparseArray, denseArray, sparseIndexes, binary::COPY );

    // both operations in one call

    HArrayUtils::buildSparseArray( sparseArray, sparseIndexes, denseArray );

