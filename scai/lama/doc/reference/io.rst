Input/Output
============

For using matrices and vectors from input files you can read the data on your own and initialize storages, matrices,
vectors with the obtained data in csr or dense format. The other way round you can get access to the stored data and
print it to file the way you like it.

Another possibility to initialize your matrices and vectors is to read data from file. 

Supported Types
---------------

 - MatrixMarket (for descriptions on the format see here__) **(matices only)**
 - SAMG format **(matrices and vectors)**
 
   - FORMATTED (ASCII)
   
   - XDR
   
   - BINARY
 
__ http://math.nist.gov/MatrixMarket/formats.html

Read from file
--------------

To initialise a matrix from file just pass the filename to the constructor or the *readFromFile* function.
For a vector pass it to the constructor.

::

   std::string matrixname = ...;
   std::string vectorname = ...;
   CSRSparseMatrix csrMatrix( matrixname );
   
   ELLSparseMatrix ellMatrix();
   ellMatrix.readFromFile( matrixname );
   
   DenseVector vec( vectorname );

Write to file
-------------

To write a matrix/vector to file call *writeToFile* with the name of the output file and the formatting. The default for just
giving a name and no formatting is binary SAMG format in internal precision for the *values* und int for *ia* and *ja*.

::

   csrMatrix.writeToFile( "matrix_out.mtx", File::MatrixMarket, File::FLOAT );
   vec.writeToFile( "vec_out", File::XDR, File::DOUBLE );
   
SAMG format
-----------

The SAMG format comes from the SAMG library of Fraunhofer SCAI and uses two files to describe a matrix or vector - one for
formatting informations, one for the data. The data can be stored in three different ways: FORMATTED, XDR, BINARY.
Formatted means the values are stored human readable in ascii, otherwise they are stored in xdr or binary format.

Matrices
^^^^^^^^

Matrix header: *.frm*
   first line:  format (f formatted, x xdr, b binary) *tab* 4 (SAMG internal version number)
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
   
