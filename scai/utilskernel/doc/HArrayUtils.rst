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
    HArrayUtils::binaryOpScalar1( result, 1.0, result, binaryOp::DIVIDE );

    // result[i] *= 2.5
    HArrayUtils::binaryOpScalar2( result, result, 2.5, binaryOp::MULT );

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
    bool ascending = true;
    HArrayUtils::sort( x, ascending );

Another possibility for sorting is usin bucket sort.

ToDo: example


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

