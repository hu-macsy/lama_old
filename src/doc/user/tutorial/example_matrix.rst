Matrix Example
==============

Multiplication of a CSRSparseMatrix and a DenseVector
-----------------------------------------------------

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
    0 & 0 & -9.3 & 4 \\
    2 & 5 & 0 & 3 \\
    2 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 2 \end{matrix}\right) *
  \left(\begin{matrix} 6 \\
    4 \\
    7 \\
    -9.3 \end{matrix}\right)
    
::

	typedef float ValueType;
	
	IndexType numRows    = 7;
	IndexType numColumns = 4;
	IndexType numValues  = 12;
	
	// VERY IMPORTANT for CSR format: ia is an offset array adding values of each row starting with 0.
	// ja stores the position of each value within a line,
	IndexType ia[] =
	{   0, 2, 3, 5, 8, 10, 10, 12};
	IndexType ja[] =
	{   0, 3, 0, 2, 3, 0, 1, 3, 0, 3, 1, 3};
	ValueType matrixValues[] =
	{   6.0f, 4.0f, 7.0f, -9.3f, 4.0f, 2.0f, 5.0f, 3.0f, 2.0f, 1.0f, 1.0f, 2.0f };
	
	// Vector values for our multiplication.
	ValueType vectorValues[] =
	{   6.0f, 4.0f, 7.0f, -9.3f};
	
	// All data has to be stored in LAMAArrays.
	const LAMAArray<IndexType> matrixIA = LAMAArray<IndexType>( numRows + 1, ia );
	const LAMAArray<IndexType> matrixJA = LAMAArray<IndexType>( numValues, ja );
	const LAMAArray<ValueType> mValues  = LAMAArray<ValueType>( numValues, matrixValues );
	const LAMAArray<ValueType> vValues  = LAMAArray<ValueType>( numColumns, vectorValues );
	
	//  Alternative code for the next 2 lines:
	//  CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>( numRows, numColumns, numValues, matrixIA, matrixJA, matrixValues );
	
	// Create a CSRStorage.
	CSRStorage<ValueType>* csrStorage = new CSRStorage<ValueType>();
	csrStorage->setCSRData( numRows, numColumns, numValues, matrixIA, matrixJA, mValues );
	
	// Allocate and fill vector for the multiplication.
	DenseVector<ValueType> vector( numColumns, 0.0 );
	vector.setValues(vValues);
	// Allocation of the result vector.
	DenseVector<ValueType> result( numRows, 0.0);
	
	// Distribution pointer are needed to construct a CSRSparseMatrix.
	lama::DistributionPtr rowDist(new lama::NoDistribution(numRows));
	lama::DistributionPtr colDist(new lama::NoDistribution(numColumns));
	
	// Allocation of the CSRSparseMatrix.
	CSRSparseMatrix<ValueType> csrMatrix(*csrStorage, rowDist, colDist);
	
	// The multiplication itself.
	result = csrMatrix * vector;
	result.writeToFile( "result.txt" , File::FORMATTED );
	
The result is:

.. math::
  \left(\begin{matrix} r_0 \\
    r_1 \\
    r_2 \\
    r_3 \\
    r_4 \\
    r_5 \\
    r_6 \end{matrix}\right) = 
  \left(\begin{matrix} -1.2 \\
    42 \\
    -102.3 \\
    4.1 \\
    2.7 \\
    0 \\
    -14.6 \end{matrix}\right)    
    
The full example program can be found here :download:`matrix.cpp <../../../examples/tutorial/matrix.cpp>`
	