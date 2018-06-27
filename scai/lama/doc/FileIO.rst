:orphan:

.. _lama_IO:

File I/O
=========

The LAMA FileIO supports the input and output of one-dimensional arrays and of (local) matrix storages into 
a file.

Note: For the IO of distributed vectors and matrices, also in mutliple files, the LAMA PartitionIO :ref:`partition_IO`
is responsible.

The FileIO supports different file types (e.g. MatrixMarket, SAMG Format, PETSc Format, usual text files) and
different file modes (binary or formatted). All value types supported by LAMA are also supported by the FileIO, i.e.
complex and real types with different precisions (float, double, long double). If the binary mode is used, array and
storage data can be saved without any loss of precision. In the formatted mode there might be a loss of precision
even if the number of significant digits can be chosen.

I/O of HArrays
--------------

Reading and writing of (heterogeneous) arrays is supported by simple read and write routines.

.. code-block:: c++

    #include <scai/lama/io/FileIO.hpp>

    using namespace scai;

    hmemo::HArray<ValueType> array;
    lama::FileIO::read( array, fileName )
    ...
    lama::FileIO::write( array, fileName )

Exceptions are thrown if files cannot be read or written. The file type is chosen by the 
suffix of the filename (e.g. ".frv", ".psc", ".txt", "mtx" ).

Type Conversions
----------------

Generally speaking, implicit type conversions are full supported. If a formatted mode is
used, there might be only the already mentioned loss of precision. 

.. code-block:: c++

    HArray<double> array1;
    FileIO::write( array1, fileName )
    HArray<float> array2;
    FileIO::read( array2, fileName )

If a binary mode is used, there might be serious problems as the file sizes are different
for different value types and the files do not contain any information about the value type
of the data stored in the file. Therefore the value type used in the file must be specified
explicitly if it is different from the value type of the array.

.. code-block:: c++

    HArray<double> array1;
    FileIO::write( array1, fileName, common::scalar::FLOAT )
    HArray<float> array2;
    FileIO::read( array2, fileName )

    HArray<double> array1;
    FileIO::write( array1, fileName )
    HArray<float> array2;
    FileIO::read( array2, fileName, common::scalar::DOUBLE )

Which solution is the preferred one, depends mainly on the application scenario in which the 
use of LAMA is embedded.

I/O of MatrixStroage
--------------------

The syntax for read and write of matrix storage is the same as for heterogeneous arrays.

.. code-block:: c++

    _MarixStorage& m = ...
    FileIO::read( m, fileName )
    ....
    FileIO::write( m, fileName )

In contrary to the heterogeneous arrays, the matrix storage classes provide also a constructor
with a single argument for the file name and provide methods for the IO.

.. code-block:: c++

    CSRStorage<double> m( fileName );

    _MarixStorage& m = ...
    m.readFromFile( fileName )
    ....
    m.writeToFile( fileName )

The rules for implicit type conversions are exactly the same as for heterogeneous arrays.

How the storage data is written into the file, is dependent on the chosen file type. Usually
a CSR ( SAMG, PETSc ) or a COO format (MatrixMarket, Text file) is used, but also the full 
dense data might be written.
The derived class from FileIO exploit type conversion routines of the matrix storage classes 
to convert between the different storage formats.

Environment Variables
---------------------

Even if the read and write routines provide addtional arguments, the file mode and the implicit
type conversions can be set by environment variables or by command line arguments.

``SCAI_IO_BINARY [bool]``

The default value is false, i.e. output is done formatted. 

``SCAI_IO_TYPE_DATA [float|double|ComplexFloat|ComplexDouble]``

This variable can be used if the data format used in the file does not match the data type in
the program. 

* Reading binary data uses a tempory buffer to read the io data of the specified type and converts it.
* Writing binary data uses a tempory buffer of the IO-type, converts values to this type before it is written.
* For formatted I/O the variable does not matter.

Be careful about the precision. Avoid conversion between complex and non-complex values as the imaginary
parts will become always zero.

You can explicitly convert a ``N x N`` complex matrix into a matrix of size ``2N x 2N`` with 
some other routines.

``SCAI_IO_TYPE_INDEX [int|long]``

In a similiar way this variable can be used to convert between 32-bit and 64-bit integer values
used for all kind of row or column indexes.

``SCAI_IO_PRECISION``

This variable can be used to set explicitly the precision of values in any formatted output.
If not set, the precision is determined by the value of the output data type.

* float, ComplexFloat: 7
* double, ComplexDouble: 12
* long, ComplexLong: 15

The precision for the aritmetic types is defined by the TypeTraits.

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

Extension for Other I/O Formats
-------------------------------

Adding support for any new file format is rather straightforward by writing
a class that derives from FileIO. Actually it should derive from CRTPFileIO
that provides already methods for the resolution from the untyped IO routine to the typed
versions, i.e. the methods with the value type as argument.

.. code-block:: c++

    class MyIO : CRTPFileIO<MyIO>, FileIO::Register<MyIO> 
    {
        static std::string createValue();   // registration value for factory
        static FileIO* create();            // create routine called for create( createValue() )

        template<typename ValueType>
        void writeStorageImpl( const MatrixStorage<ValueType>& storage, const std::string& fileName );

        template<typename ValueType>
        void readStorageImpl( MatrixStorage<ValueType>& storage, const std::string& fileName );

        template<typename ValueType>
        void writeArrayImpl( const hmemo::HArray<ValueType>& array, const std::string& fileName );
        __attribute( ( noinline ) );

        template<typename ValueType>
        void readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName );
    }

The typical use of such an IO File handler class would be as follows:

.. code-block:: c++

    void function( ..., _MatrixStorage& storage, ... )

       MyIO myIO;
       myIO.readStorage( storage, fileName );
       // e.g. calls myIO->readStorageImpl<double>( any, fileName ), if storage->getValueType() == DOUBLE

Any IO-Handler is intended to register itself in the FileIO factory. For the registration
the file suffix is used as key to create the corresponding IO handler. 
The following code shows a typical example in the LAMA core where 
the factory is used to call a virtual routine that results in calling a corresponding
method of the IO handler.

.. code-block:: c++

    void function( ..., const std::string& inFileName, _MatrixStorage& storage )
    {
        ...
        std::string suffix = FileIO::getSuffix( inFileName );

        if ( FileIO::canCreate( suffix ) )
        {
            // okay, we can create derived FileIO object by factory
    
            common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );
    
            SCAI_LOG_INFO( "Got from factory: " << *fileIO  )

            fileIO->any_virtual_fn( ..., storage, ... )
        
         }
         ...
     }


