Matrix Example
--------------

Multiplication of a CSRSparseMatrix and a DenseVector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we simply multiply a sparse matrix in CSR format and a dense vector.
Our multiplication looks like this:   

.. math::

  \left(\begin{matrix} r_0 \\
    r_1 \\
    r_2 \\
    r_3 \\
    r_4 \\
    r_5 \\
    r_6 \end{matrix}\right) =
  \left(\begin{matrix} 6 & 0  & 0 & 4 \\
    7 & 0 & 0 & 0 \\
    0 & 0 & -9 & 4 \\
    2 & 5 & 0 & 3 \\
    2 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 2 \end{matrix}\right) *
  \left(\begin{matrix} 6 \\
    4 \\
    7 \\
    -9 \end{matrix}\right)

.. literalinclude:: ../../lama/examples/tutorial/matrix.cpp 
   :language: c++
   :lines: 35-122
    
The result is:

.. math::

  \left(\begin{matrix} r_0 \\
    r_1 \\
    r_2 \\
    r_3 \\
    r_4 \\
    r_5 \\
    r_6 \end{matrix}\right) = 
  \left(\begin{matrix} 0 \\
    42 \\
    -99 \\
    5 \\
    3 \\
    0 \\
    -14 \end{matrix}\right)    
    
The full example program can be found here :download:`matrix.cpp <../../lama/examples/tutorial/matrix.cpp>`
	
Optimization of Setting Up the Matrix and Vector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the above example program the matrix data was built in several steps:

- copying the raw data to heterogeneous arrays (to allow operations on different devices)
- copying the heterogenoeous arrays to an CSR storage object (storage is compounded of arrays for sparse pattern and non-zero values)
- copying the CSR storage object into a CSR matrix object (to allow the distribution of matrix data and parallel operations).

.. code-block:: c++

    HArray<IndexType> csrIA( numRows + 1, rawIA );
    HArray<IndexType> csrJA( numValues, rawJA );
    HArray<ValueType> csrValues( numValues, rawValues );
    CSRStorage<ValueType> csrStorage ( numRows, numColumns, csrIA, csrJA, csrValues );
    CSRSparseMatrix<ValueType> csrMatrix( csrStorage );

One copy, here the first one, is mandatory as otherwise data cannot be modified and managed on its own.
The two other copies are not really needed if we do not care about the data after it has been copied.
By using the move semantics of C++11, it is possible to avoid these copy steps as follows:

.. code-block:: c++

    CSRStorage<ValueType> csrStorage ( numRows, numColumns, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );
    CSRSparseMatrix<ValueType> csrMatrix( std::move( csrStorage ) );

Please note that the move operations leave the heterogeneous arrays ``csrIA``, ``csrJA``, and ``csrValues``
as well as the CSR storage ``csrStorage`` as empty containers.

By comparing the pointers of the allocated data it is possible to verify the correctness of the move semantic:

.. code-block:: c++

   const ValueType* ptr1 = ReadAccess<ValueType>( csrValues ).get();
   ...
   const ValueType* ptr2 = ReadAccess<ValueType>( csrMatrix.getLocalStorage().getValues() ).get();

   SCAI_ASSERT_EQ_ERROR( ptr1, ptr2, "data has not been moved correctly" )

The full example program can be found here :download:`matrix_opt.cpp <../../lama/examples/tutorial/matrix_opt.cpp>`

Please keep in mind that the move semantic is only availabe for data of same value type and same storage format.

