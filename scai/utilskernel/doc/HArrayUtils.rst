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

Reduction Operations
--------------------

.. code-block:: c++

    HArray<double> x;

    double sum = scan( x );
    double sum = scan( x );

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

Sparse Arrays
-------------

The following operations can be used to build a sparse array from a dense array
and vice versa.

.. code-block:: c++

    HArray<double> denseArray;
    HArray<IndexType> sparseIndexes;

    HArrayUtils::buildSparseIndexes( sparseIndexes, denseArray );

