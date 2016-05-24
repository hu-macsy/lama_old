.. _lama_Matrix:

Matrix
======

The class ``Matrix`` is a generic mathematical matrix. The index and data arrays are internally stored in a ``HArray`` out of :ref:`hmemo <scaihmemo:main-page_hmemo>` so a ``Matrix`` can transparently used on every device.

A specific representation of a matrix can be a DenseMatrix or a sparse matrix of a specific storage type. A DenseMatrix saves all entries of the m times n matrix, while a sparse matrix only stores the non-zero entries in a well defined format. LAMA actually preserves the following sparse matrix formats (for detailed descriptions on the formats refer to the subpages of :ref:`sparsekernel <scaisparsekernel:main-page_sparsekernel>`:

 - COO (Coordinate)
 - CSR (Compressed Sparse Row)
 - DIA (Diagonal)
 - ELL (ELLPACK-R)
 - JDS (Jagged Diagonal Storage)

Constructors
------------

The class ``Matrix`` is an abstract class that can be used for generic algorithm formulation.
For instantiating a matrix variable you need to call the constructor one of the following templated matrix classes, that are a specific representation of a matrix holding the matrix entries:

 * COOSparseMatrix
 * CSRSparseMatrix
 * DIASparseMatrix
 * DenseMatrix
 * ELLSparseMatrix
 * JDSSparseMatrix

For creating a new Matrix you need these major things:
 * the number of rows of the matrix
 * the number of columns of the matrix
 * the data of the matrix either in dense or csr representation

For distributed matrices the number of rows/columns can be substituted by a ``Distribution`` (holding the number of rows/columns as row/column distribution strategy). For defining a Distribution, please refer to :ref:`this <scaidmemo:main-page_dmemo>` page.

.. code-block:: c++

  // for later use:
  int numRows    = 6;
  int numColumns = 4;
  dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr( "MPI" ) );
  dmemo::DistributionPtr rDist( dmemo::Distribution::getDistribution( "BLOCK", comm, numRows, 1.0 ) );
  dmemo::DistributionPtr cDist( dmemo::Distribution::getDistribution( "BLOCK", comm, numColumns, 1.0 ) );
  common::ContextPtr cudaContextPtr = common::Context::getContextPtr( common::context::CUDA );

  // empty (not initialized) float matrices (with context, distribution, or both)
  DenseMatrix<float> a();                      // without size
  DenseMatrix<float> b( numRows, numColumns ); // local matrix with size 6x4
  DenseMatrix<float> c( rDist, cDist );        // block distributed matrix with global size 6x4
  DenseMatrix<float> square( rDist );          // block distributed matrix with global size 6x6

  // initialize with csr representation
  int nnz = 13
  int ia[] = {0, 2, 4, 6, 9, 11, 11, 13};
  int ja[] = {0, 3, 1, 0, 2, 3, 3, 0, 1, 0, 3, 1, ,3};
  double values[] = {6.0, 4.0, 0.0, 7.0, 9.0, 4.0, 3.0, 2.0, 5.0, 2.0, 1.0, 1.0, 2.0};
  DenseMatrix<double> csrFormat( numRows, numColumns, nnz, &ia, &ja, &values );


  // copy constructor (also works with general Matrix 'a')
  DenseMatrix<double> aCopy   ( a );
  DenseMatrix<double> aRedist ( a, rDist, cDist ); // a with new Distributions



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
:doc:`IO`. All data read in is stored in CSRStorage and converted to the target format afterwards.

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
       