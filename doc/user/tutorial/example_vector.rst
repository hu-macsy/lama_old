Scalars and Vectors
===================

Multiplication of a Scalar and a DenseVector
--------------------------------------------

In this example we simply multiply a scalar with a dense vector.  

::

    typedef double ValueType;     

    Scalar singleValue( 2.0 );

    const ValueType inputData[] = { 1.0, 2.0, 3.0, 4.0 };

    DenseVector<ValueType> sequenceOfValues( 4, inputData );

    sequenceOfValues = singleValue * sequenceOfValues;

    sequenceOfValues.writeToFile( "vector" , File::FORMATTED );

The full example program can be found here :download:`vector.cpp <../../../examples/tutorial/vector.cpp>`
