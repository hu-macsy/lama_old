.. _lama_SetStorage:

Setting a Storage
=================

For writing algorithms you can use the generic data structur ``Matrix`` without the decision about the internal storage format, but when implementing dedicated applications you need to choose one. Choosing the right storage format will influence the execution performance and should be driven by two main factors:

1. the structure of the matrix
2. the exection context (device)

There are several factors that need to be taken into account for choosing the right format. Generally it can be seperated in the following questions?

* Does it have a lot of zero values? 
	* YES: use one of the sparse matrices
	* NO:  use the Dense format
* If sparse: 
	* If your matrix is strictly banded: use the DIA format
	* If your matrix has nearly the same number of nonzero entries per row: use the ELL format
	  * If there are few rows with a lot more nonzeros: use the JDS format (ELL with sorted rows)
	* If your matrix has a irregular pattern: use the CSR format
	  * If you already use random access to the matrix elements: use the COO format

For detailed information about the storage formats, please refer to :ref:`scaisparsekernel:main-page_sparsekernel`.

You can set or use a specific storage format in multiple ways, e.g if you want a CSR storage format storing the matrix entries as float:

1. By using a specific matrix type with a specific data type (result: you have a matrix variable):

.. code-block:: c++

	CSRSparseMatrix<float> csrMatrix( ... ); // calling a constructor, empty: default constructor

2. By calling the matrix factory with a create key describing matrix and data type (result: you have a matrix pointer variable):

.. code-block:: c++
	
	MatrixCreateKeyType matrixType( Format::CSR, common::scalar::FLOAT );
	lama::MatrixPtr matrixPtr( Matrix::create( MatrixCreateKeyType ( matrixType ) ) ); // calling the default constructor
