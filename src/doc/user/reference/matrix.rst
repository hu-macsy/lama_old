Matrix
======

The class Matrix is a generic mathematical matrix. A specific representation of a matrix can be a DenseMatrix or a sparse
matrix of a specific storage type (e.g. compressed sparse rows (CSR): CSRSparseMatrix). A DenseMatrix saves all entries
of the m times n matrix, while a sparse matrix only stores the non-zero entries in a well defined format.

storages
--------

LAMA preserves the following sparse matrix formats:

 - compressed sparse row (CSR)
 - ELLpack (ELL)
 - jagged diagonal storage (JDS)
 - coordinate (COO)
 - diagonal (DIA)
 
For a detailed description on the different storage types, see :doc:`storage`.

You can create the specific matrix types by the same constructor giving the number of rows and columns:

.. code-block:: c++

   int m = 7, n = 4;
   CSRSparseMatrix<ValueType> mCSR ( m, n );
   ELLSparseMatrix<ValueType> mELL ( m, n );
   ...
   DenseMatrix<ValueType> mDense ( m, n );
   
Or you create it from a matrix storage. As example we fill a CSRStorage with CSR data to create a CSRSparseMatrix:

.. code-block:: c++

    IndexType ia[] = { 0, 2, 3, 5, 8, 10, 10, 12 };
    IndexType ja[] = { 0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3 };
    ValueType values[] = { 6.0f, 4.0f, 7.0f, -9.3f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };
    int numValues = 12;
    
    const LAMAArray<IndexType> mIA = LAMAArray<IndexType>( m + 1, ia );
    const LAMAArray<IndexType> mJA = LAMAArray<IndexType>( numValues, ja );
    const LAMAArray<ValueType> mValues  = LAMAArray<ValueType>( numValues, values );
    
    CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( m, n, numValues, mIA, mJA, mValues );
    
    CSRSparseMatrix<ValueType> csrMatrix( *csrStorage );
    
Conversion between the different matrix formats is preserved internally, so you can create each matrix type out of
another:

.. code-block:: c++

   ELLSparseMatrix<ValueType> ellMatrix( csrMatrix );
   
Of course you can create a matrix from file. For a detailed description on supported file formats please refer to
:doc:`io`. All data read in is stored in CSRStorage and converted to the target format afterwards.

.. code-block:: c++

   JDSSparseMatrix<ValueType> jdsMatrix( "matrix.mtx" );
   
As with vectors you can perform matrix expressions as addition, substraction, multiplication with other matrices or
vectors and scaling with a scalar in text-book syntax:

.. code-block:: c++

   Scalar s( 2.0 );
   DenseVector x( m, 1.0 );
   DenseVector y( m, 3.0 );
   
   ellMatrix = s * ellMatrix;
   mELL = mELL + ellMatrix;
   
   DenseVector r = csrMatrix * x;
   r = csrMatrix * x + y;
   r = csrMatrix * x + 2.0 * y;
   r = csrMatrix * x + y * 2.0;
   
Useful functions on a matrix to get to know sizes are:

.. code-block:: c++

   m = jdsMatrix.getNumRows();
   n = jdsMatrix.getNumColumns();
   numValues = jdsMatrix.getNumValues();
   
You also can get the matix maximum norm by:

.. code-block:: c++

   s = csrMatrix.maxNorm();
   
To save a matrix to file use 'writeToFile':   
   
.. code-block:: c++
 
   ellMatrix.writeToFile( "output.mtx", File::MatrixMarket, File::FLOAT );
       