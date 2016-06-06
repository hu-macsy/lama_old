.. _lama_SetValueType:

Setting/Using a ValueType
=========================

For writing algorithms you can use the generic data structures ``Scalar``, ``Vector`` and ``Matrix`` without the decision about a datatype of your vector or matrix, but when implementing dedicated applications you need to choose one.

``Scalar`` can hold every data type without a special usage. It saves the data type it is given. So the user has nothing special to do. Used in expressions it adaptes the data type of the used matrices or vectors. Initilizing a Scalar is handled the same as a the basic arithmetic types, or by explicit assignment - complex data type always need the explicit assigment.

* int    for single digital value, e.g. 2
* float  for floating point value with f, e.g. 2.53f
* double for floating point value, e.g. 2.53
* ComplexFloat,  e.g. ComplexFloat (3.2, 1.4)
* ComplexDouble, e.g. ComplexDouble(3.2, 1.4)

``Vector`` and ``Matrix`` are abstracted classes having no data type. Only their inherited classes ``DenseVector``, ``CSR-/ELL-/...SparseMatrix`` and ``DenseMatrix`` have a specific data type you can define as template type in rectangular brackets (``<>``), e.g. ``DenseVector<float>``, ``CSRSparseMatrix<double>`` or ``DenseMatrix<ComplexFloat>``.
