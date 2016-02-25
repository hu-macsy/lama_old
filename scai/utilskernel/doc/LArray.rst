The Local Array LArray
======================

The array class **LArray** stands for a local array and is derived from the class *HArray*. Regarding
the member variables it is exactly the same, but it offers many addtional operators.

.. code-block:: c++

    template<typename T>
    class LArray : HArray<T>
    { 
    }

For the implementation of the operators kernel functions implemented on different devices
are used.

.. code-block:: c++

  LArray<double> A;
  input( A );
  LArray<float> B( A );   // implicit conversion
  ...
  double x = A[0];
  A[1] = x * 2.0;
  A[2] = A[1];

  LArray<IndexType> IND;
  A[IND] = C;  !
  C = B[IND];

  A = 2.0 * B;
  A = A * B;   // componentwise multiplication
