.. _file_format:

Supported File Formats
======================

The following table lists all supported file formats by the current LAMA version.
All file formats allow to store at least a dense vector or a sparse matrix
(either CSR or COO format). And for each supported format it is possible to use
the single or independend mode for distributed data.

.. list-table:: Supported File Formats
   :widths: 10 10 10 10 10 10 10
   :header-rows: 1

   * - Name
     - Lama
     - MatrixMarket
     - Text
     - SAMG
     - PETSc
     - Matlab
   * - Suffixes
     - lmf
     - mtx
     - txt
     - frm/frv
     - psc
     - mat
   * - Collective Mode
     - yes
     - no
     - no
     - no
     - no
     - no
   * - File Mode
     - Binary
     - Text
     - Text
     - Binary/Text
     - Binary
     - Binary
   * - Multiple items
     - yes
     - no
     - no
     - no
     - yes
     - yes
   * - sparse vector
     - yes
     - yes
     - no
     - no
     - no
     - no
   * - dense matrix
     - yes
     - yes
     - no
     - no
     - no
     - no
   * - grid vector
     - yes
     - no
     - no
     - no
     - no
     - yes

Here are some general remarks:

 - If a format does not support a sparse vector, it will be written as a 
   dense vector that might require some more file space.
 - If a format does not support grid vector, it is just written as a dense vector.
   Number of dimensions and sizes of each dimension is not available.
 - If a format does not support a dense matrix it is written as a  sparse matrix.
 - If a format does not support the collective mode, LAMA uses automatically the serial 
   mode for read/write operations.
 - Matrix-Market, Text and SAMG format only allow to write one item (vector or matrix) into 
   a file. An append mode is not supported for the corresponding file objects.

Supported File Types
--------------------

The decision about the file type is taken by the suffix of the file name:
Currently, the following file types are supported

 - MatrixMarket (for description on the format see |MM|), for suffix ".mtx"

 - SAMG format (see below), for suffix ".frm" (matrix) or ".frv" (vector)
 
   - FORMATTED (ASCII)
   
   - BINARY

 - PETSC format (binary format), for suffix ".psc"

 - Text format (pure ASCII format), for suffix ".txt"

 - Level 5 MAT-File format of Matlab, for suffix ".mat" 

.. |MM| raw:: html

   <a href="http://math.nist.gov/MatrixMarket/formats.html" target="_blank"> here </a>

Conversion from one file type to another file type is rather simple, just read the matrix/vector from one file
and write it with its new extension to another file.

.. code-block:: c++

    _MatrixStorage& m = ...
    m.readFromFile( "matrix_3D27P_100_100_100.txt" )
    m.writeToFile( "matrix_3D27P_100_100_100.mtx" )

Here are some remarks:

 * The matrix type, e.g. CSR, DIA, ELL, JDS, does not matter when reading or writing matrix data.
   There is always an implicit conversion when reading or writing the data. Nevertheless the CSR format
   is preferred as it has usually the minimal overhead.
 * The value type, e.g. float, double, ComplexFloat, ComplexDouble is taken over if the binary mode is used,
   i.e. there is no loss of precision. In the formatted output, the number of significant digits depends on
   the value type, but there may be a certain loss of precision. Implicit type conversion is supported but
   should be used rather carefully.
 * Usually a certain file type supports the formatted or the binary mode. Only the SAMG format supports both modes.
 * Some formats do not store for a matrix storage the number of columns explicitly. Here the number of columns
   is determined by the maximal column index that appears in the data.

SAMG format
-----------

The SAMG format comes from the |SAMG| library of Fraunhofer SCAI and uses two files to describe a matrix or vector - 
one header file with general information ( mode, size), one for the data. 
The data can be saved in both modes, either BINARY or FORMATTED.

.. |SAMG| raw:: html

   <a href="https://www.scai.fraunhofer.de/de/geschaeftsfelder/schnelle-loeser/produkte/samg.html" target="_blank"> SAMG </a>

Matrices
^^^^^^^^

Matrix header: *.frm*
   first line:  mode (f formatted, b binary) *tab* 4 (SAMG internal version number)
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
   first line: mode (f formatted, x xdr, b binary)
   second line: number of values (nv)
   third line: size of value type (in most cases: 4 for float, 8 for double)
   
Vector data: *.vec*
   nv lines with values (one value per line)

Level 5 MAT-File Format
-----------------------

This file format can be used to exchange data between LAMA and MATLAB applications.

The following examples shows a MATLAB code that generates a random matrix that is written to a file.

.. code-block:: c++

   >> mat = sprand( 6, 4, 0.3 )

   mat =

      (2,1)       0.7952
      (3,1)       0.4898
      (1,2)       0.3816
      (5,2)       0.6463
      (2,3)       0.1869
      (1,4)       0.7655
      (4,4)       0.4456

   >> save /home/brandes/MAT/sp6x4.mat mat

This matrix can be read in LAMA as follows:

.. code-block:: c++

   scai::lama::CSRSparseMatrix<double> mat( "/home/brandes/MAT/sp6x4.mat" );
   mat = 2 * mat;
   mat.writeToFile( "/home/brandes/MAT/sq6x4.mat" );

.. code-block:: c++

   >> load /home/brandes/MAT/sq6x4.mat
   >> LAMA

   LAMA =

      (2,1)       1.5904
      (3,1)       0.9795
      (1,2)       0.7631
      (5,2)       1.2926
      (2,3)       0.3737
      (4,4)       0.8912
      (1,4)       1.5310

When using this file format the following issues should be considered:

 - LAMA can only read and write one single data element from a MATLAB file.
 - The name of the data element is is ignored when reading the element and 
   each written element gets the name "LAMA".
 - As the data type is stored for each element in the file, the SCAI_IO_TYPE
   is ignored, i.e. each array/storage is written exactly in the format it is
   and there might be an implicit type conversion during the read.
 - Only the SCAI_IO_TYPE=PATTERN will be handled and in this case no sparse matrix
   values are written, only the row and column indexes
 - The data types LongDouble and ComplexLongDouble are not really supported
   by the Level 5 MAT-File format but are full supported here by using a reserved
   value of the MAT-File Data Types.

Text Format
-----------

When using the text format a (sparse) matrix is saved in the COO format where
each line contains the row index, column index and value of the non-zero entry.

.. code-block:: c++

      1 0  1.5904
      2 0  0.9795
      0 1  0.7631
      4 1  1.2926
      1 2  0.3737
      3 3  0.8912
      0 3  1.5310

When using this file format the following issues should be considered:

 - The number of non-zero entries is given by the number of lines of the file
 - The number of rows and columns of the matrix is not stored explicitly but is
   determined by the maximal row and column index when reading the file.

A (dense) vector is saved in a text file by one entry for each value of the vector.

.. code-block:: c++

     0.51 
     0.43
     0.31 

Using the text format is not recommended for efficient I/O but might be very useful
for testing and developping. Furthermore, it might be a convenient way to exchange
data with other applications. Here is an example of how to use this format to
exchange data with MATLAB applications.

.. code-block:: matlab

 [i,j,val] = find( matrix )          [i,j,val] = find(matrix)
 data_dump = [i, j, val] A           fid = fopen( 'data.txt', 'w' )
 save -ascii data.txt data_dump      fprintf( fid, "%d %d %f\n", [i,j,val] )
                                     flose( fid )

.. code-block:: matlab

  data_dump = importdata( 'data.txt' )      load data.txt
  matrix = spconvert( data_dump )           matrix = spconvert( data )

LAMA Format
-----------

This is a proprietary binary format and offers in a certain sense
the advantages of all other formats.

- uses CSR format with sizes instead offsets
- stores information about the data type 
- supports sparse vector
- supports grid data

The only disadvantages:

- not readable
- not supported by other applicatins

do not really matter as conversions between all other formats are supported.

