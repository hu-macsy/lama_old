:orphan:

.. _lama_IO:

Input/Output
============

For creating matrices and vectors from input files you can read the data from file in supported format types.

If you have your own special file format you have to read the data on your own and initialize storages, matrices, vectors with the obtained data in csr or dense format. The other way round you can get access to the stored data and print it to file the way you like it.

Supported Types
---------------

 - MatrixMarket (for description on the format see |MM|)

 - SAMG format (see below)
 
   - FORMATTED (ASCII)
   
   - BINARY

   - XDR
 
.. |MM| raw:: html

   <a href="http://math.nist.gov/MatrixMarket/formats.html" target="_blank"> here </a>

SAMG format
-----------

The SAMG format comes from the |SAMG| library of Fraunhofer SCAI and uses two files to describe a matrix or vector - one for formatting informations, one for the data. The data can be stored in three different ways: FORMATTED, BINARY, XDR. Formatted means the values are stored human readable in ASCII, otherwise they are stored in binary or xdr format.

.. |SAMG| raw:: html

   <a href="https://www.scai.fraunhofer.de/de/geschaeftsfelder/schnelle-loeser/produkte/samg.html" target="_blank"> SAMG </a>

Matrices
^^^^^^^^

Matrix header: *.frm*
   first line:  format (f formatted, b binary, x xdr) *tab* 4 (SAMG internal version number)
   second line: *tab tab* number of values (nv) *tab* number of rows (nr) *tab* 22 (SAMG internal: symmetry information) *tab* 1 (SAMG internal: number of unknowns ) *tab* 0 (SAMG internal)   

.. 22: unsymmetric, not equal sums of row

Matrix data: *.amg*
   one value per line:
   nr lines with ia data
   nv lines with ja data
   nv lines with values
   
Vectors
^^^^^^^

Vector header: *.frv*
   first line: format (f formatted, x xdr, b binary)
   second line: number of values (nv)
   third line: size of value type (in most cases: 4 for float, 8 for double)
   
Vector data: *.vec*
   nv lines with values (one value per line)

Read from file
--------------

To initialise a ``Matrix`` or ``Vector`` from file just pass the filename to the constructor or the *readFromFile* function.

.. code-block:: c++

   std::string matrixname = ...;
   std::string vectorname = ...;
   CSRSparseMatrix<float> csrMatrix( matrixname );
   
   ELLSparseMatrix<double> ellMatrix();
   ellMatrix.readFromFile( matrixname );
   
   DenseVector<float> vec1( vectorname );

   DenseVector<double> vec2();
   vec2.readFromFile( vectorname );

Write to file
-------------

To write a ``Matrix`` or ``Vector`` to file call *writeToFile* with the name of the output file and the formatting. The default for just giving a name and no formatting is binary SAMG format in internal precision for the *values* und int for *ia* and *ja*.

.. code-block:: c++

   csrMatrix.writeToFile( "matrix_out.mtx", File::MATRIX_MARKET, File::FLOAT );
   vec.writeToFile( "vec_out.frv", File::SAMG_FORMAT, File::DOUBLE, /*binary*/ true ); // binary SAMG format
   
Possible file formats are ``File::SAMG_FORMAT`` and ``File::MATRIX_MARKET``.

Possible data types are ``common::scalar::INDEX_TYPE`` (int), ``common::scalar::FLOAT``, ``common::scalar::DOUBLE``, ``common::scalar::COMPLEX``(ComplexFloat), ``common::scalar::DOUBLE_COMPLEX``, ``common::scalar::LONG_DOUBLE_COMPLEX`` or ``common::scalar::INTERNAL`` for the internal representation of the data.
