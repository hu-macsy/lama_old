Vector
======

The class Vector is a generic mathematical vector. A specific representation of a vector is a DenseVector, holding the
full vector entries.

You can create a DenseVector by passing the number of elements (size) and a sequence of input values to the constructor,
by passing the size and one value assigned to each entrie of the vector or by reading it from file. 

.. code-block:: c++

    const ValueType inputData[] = { 1.0, 2.0, 3.0, 4.0 };
    DenseVector<ValueType> x( 4, inputData );

    DenseVector<ValueType> y( "in.txt" );

    DenseVector<ValueType> z( 4, 0.0 );

Having vectors and scalars (as Scalar or value) you can perform vector addition, subtraction and scaling with a scalar
in text-book syntax:

.. code-block:: c++

    Scalar s( 2.0 );
    x = s * x;
    
    z = x + y;
    z = x * 2.0 + y;
    z = 2.0 * x + y;
    z = x + 1.0 * y;
    z = x + y * 1.0;
    
    z = y * 2.0;
    z = y / 2.0;
    
and also combined with the assign operator:

.. code-block:: c++

    z += x;
    z += 2.0 * x;
    z += x * 2.0;

    z -= x;
    z -= 2.0 * x;
    z -= x * 2.0;
    z *= 3.0;
    z /= 1.5;

You can assign one value to the whole vector also by '=' (size is kept):

.. code-block:: c++

    x = 1.0;
    y = 2.0;
    
The dot product of two vectors is expressed as function (so there is no misunderstanding):

.. code-block:: c++

   s = x.dotProduct( y );

You can get the L1-, L2-, maximum-norm of an vector by:
   
.. code-block:: c++ 
   
   s = x.l1Norm();
   s = x.l2Norm();
   s = x.maxNorm();

Other useful functions on a vector are:

.. code-block:: c++ 

   IndexType length = x.size(); // getting the size of a vector
   
   Scalar maximum = x.max(); // getting the maximum value of all entries
   Scalar minimum = x.min(); // getting the minimum value of all entries
   
   // writing a vector to file as formatted output in double precision
   z.writeToFile( "vector.txt", File::FORMATTED, File::DOUBLE);