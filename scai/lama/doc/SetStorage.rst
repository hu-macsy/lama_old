.. _lama_SetStorage:

Setting a Storage
=================

For writing algorithms you can use the generic data structur ``Matrix`` without the decision about the internal storage format, but when implementing dedicated applications you need to choose one. Choosing the right storage format will influence the execution performance and should be driven by two main factors:

1. the structure of the matrix
2. the exection context (device)

The structure of the matrix is the first factor that needs to be taken into account. It can be seperated in the following questions?

* Does it have a lot of zero values? 
	YES: Use one of the sparse matrices.
	NO: Use ``DenseMatrix``.
* If sparse: 

