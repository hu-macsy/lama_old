:orphan:

.. _lama_IO:

Input/Output
============

For creating matrices and vectors from input files you can read the data from file in supported file types.

If you have your own special file typee you have to read the data on your own and initialize storages, matrices, 
vectors with the obtained data in csr or dense format. The other way round you can get access to the stored data 
and print it to file the way you like it.

I/O of HArrays
==============

.. code-block:: c++

    HArray<ValueType> array;
    FileIO::read( array, fileName )
    ...
    FileIO::write( array, fileName )

I/O of MatrixStroage
====================

.. code-block:: c++

    _MarixStorage& m = ...
    FileIO::read( m, fileName )
    ....
    FileIO::write( m, fileName )

.. code-block:: c++

    _MarixStorage& m = ...
    m.readFromFile( fileName )
    ....
    m.writeToFile( ... )

.. code-block:: c++

    CSRStorage<double> m( fileName );

I/O of Vectors
==============

In contrary to a single array, the vector data might be distributed among the different processors.

.. code-block:: c++

    Vector& v = ...
    v.readFromFile( fileName )
    ....
    v.writeToFile( ... )

When writing a distributed vector to a single file, data is gathered on the master process that writes the complete
data into the file. 

When reading a vector from a single file, data is only read by the master process. The distribution set for the
vector specifies exactly this mapping. Nevertheless, the data might be redistributed later within the application as
required.

A vector can be written to and read from mutliple file, one for each processor. This is exactly the case when
the fileName contains the substring "%r" that is replaced with "<rank>.<size>" of the communicator.

.. code-block:: c++

    Vector& v = ...
    ....
    v.writeToFile( "data_%r.mtx" )

.. code-block:: c++

    mpirun -np 2 <appl>  -> write files data_0.2.mtx and data_1.2.mtx
    mpirun -np 3 <appl>  -> write files data_0.3.mtx, data_1.3.mtx and data_2.3.mtx

But be careful. If the distribution of the vector is not a block distribution, the mapping information
is lost. When reading the vector from multiple files, the distribution of the vector will be set
implicitly to a corresponding general block distribution.

I/O of Distributions
====================

If a vector or a matrix is stored in multiple files, the information about the mapping is lost 
if the distribution is not a (general) block distribution.

Therefore the mapping itself can be written to a single or mutliple files. 

A single file contains for each entry the owner.

.. code-block:: c++

   0   0   1    1   0   0   2   2   3   3   2   2   3   3

A multiple file contains for each processor the owned global indexes.

.. code-block:: c++

   0 :   0  1   4   5
   1 :   2  3   6   7
   2 :   8  9  12  13
   3 :  10 11  14  15

.. code-block:: c++

   PartitionIO::write( distribution, "ownwers.mtx" )

I/O of Matrices
===============

Writing and reading a matrix to a single file is done in the same way as for a vector.

.. code-block:: c++

    m.writeToFile( "matrix.mtx" )  -> all data is gathered on the master process and written
    m.readFromFile( "matrix.mtx" )  -> all data is on the master process

    m.writeToFile( "matrix_%r.mtx" ) -> each processor reads a local part of the matrix
    m.readFromFile( "matrix_%r.mtx" ) -> each processor reads a local part of the matrix

Reading a matrix 

.. code-block:: c++

    DistributionPtr dist = ...
    m.readFromFile( "matrix_%r.mtx", dist )

.. code-block:: c++

    m.readFromFile( "matrix_%r.mtx", "owners.mtx" )
    m.readFromFile( "matrix_%r.mtx", "myIndexes%r.mtx" )

In contrary to a vector, a partitioned matrix might still contain the info about its
distribution. This is the case if the data is stored in a sparse format and the first
column index of a row is the diagonal element. As column indexes are still global, the
array of first column indexes for each row is the same as the global indexes of each partition
stored in a partitioned mapping file.

.. code-block:: c++

    m.readFromFile( "matrix_%r.mtx", "" )

Consider the following example of a 16 x 16 matrix:

.. code-block:: c++

   matrix        owner, local index   

   0 0 4             0  0          
   1 1 4             0  1 
   2 2 4             1  0
   3 3 4             1  1
   4 4 4             0  2
   5 5 4             0  3
   6 6 4             1  2
   7 7 4             1  3
   8 8 4             2  0
   9 9 4             2  1
   10 10 4           3  0 
   11 11 4           3  1
   12 12 4           2  2
   13 13 4           2  3
   14 14 4           3  2
   15 15 4           3  3
   0 1 -1
   ....

.. code-block:: c++

   matrix_1.0.txt   matrix_1.4.txt    matrix_2.4.txt    matrix_3.4.txt

   0 0 4            0  2  4            0  8  4            0  10  4
   1 1 4            1  3  4            1  9  4            1  11  4
   2 4 4            2  6  4            2 12  4            2  14  4
   3 5 4            3  7  4            3 13  4            3  15  4
   0 1 -1           0  1 -1            0  4  -1           0  6  -1
   ....             ...                ...                ...

Supported File Types
--------------------

The decision about the file type is taken by the suffix of the file name:
Currently, the following file types are supported

 - MatrixMarket (for description on the format see |MM|), for suffix ".mtx"

 - SAMG format (see below), for suffix ".frm" (matrix) or ".frv" (vector)
 
   - FORMATTED (ASCII)
   
   - BINARY

 - PETSC format (binary format), for suffix ".psc"

 - MATLAB format (pure ASCII format), for suffix ".txt"

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

Environment Variables
---------------------

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

Example Program
---------------

In the direcory ``eamples/io`` two example programs are provided that allow nearly each kind of conversion between
different file formats: one for matrices and one for vectors.

.. code-block:: bash

   matrixConvert <infile_name> <outfile_name> [options]
   vectorConvert <infile_name> <outfile_name> [options]

   SCAI_TYPE=[float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLongDouble

   SCAI_IO_BINARY=flag[:bool]
   SCAI_IO_TYPE_DATA=string[:float|double|ComplexFloat|ComplexDouble]
   SCAI_IO_TYPE_INDEX=[int]
   SCAI_IO_PRECISION=[n:int]
   SCAI_IO_APPEND=flag

Here are some examples:

.. code-block:: bash

   matrixConvert mhd1280b.mtx mhd1280b.frm --SCAI_TYPE=ComplexDouble

This converts a complex matrix (MatrixMarket) to the binary SAMG format.

.. code-block:: bash

   matrixConvert Emily_923.mtx Emily_923.psc 

This converts a double matrix (MatrixMarket) to the binary PETSC format.

.. code-block:: bash

   matrixConvert matrix.frm file.psc 
   vectorConvert matrix.frv file.psc --SCAI_IO_APPEND=True

This converts a double matrix and a double vector into one single binary PETSC file.

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


