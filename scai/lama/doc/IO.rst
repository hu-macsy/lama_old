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


