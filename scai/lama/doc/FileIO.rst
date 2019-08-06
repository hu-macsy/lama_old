.. _lama_IO:

FileIO
======

FileIO is an abstract class, all other IO classes derive from it.

FileIO Class Hierarchy
----------------------

.. figure:: _images/single_io.*
    :width: 600px
    :align: center
  
    Single I/O of Distributed Data.

Each derived class must implement at least these pure (serial) methods:

 * open and close of a file
 * write an array (dense vector) into the file
 * read an array (dense vector) from the file
 * write a coordinate (sparse) matrix into the file
 * read a coordinate (sparse) matrix into the file

.. code-block:: c++

    #include <scai/lama/io/MatrixMarketIO.hpp

    using namespace scai;

    hmemo::HArray<double> values = ...
    CSRStorage<double> csr = ...
    
    MatrixMarketIO file;

    // output of values, csr             // input of values, csr

    file.open( "vector.mtx", "w" );      file.open( "vector.mtx", "r" );
    file.writeArray( values );           file.readArray( values );
    file.close();                        file.close();
    file.open( "matrix.mtx", "w" );      file.open( "matrix.mtx", "r" );
    file.writeStorage( csr );            file.readStorage( csr );
    file.close();                        file.close();

The above code is exactly the same for all other supported file types.

FileIO Factory
--------------

As seen before, most I/O operations can be written in such a way that they
are independent of the underlying code. 

.. code-block:: c++

    FileIO& file;

All I/O classes can register themselves into the FileIO-Factory. This allows to 
create dynamically a certain FileIO object of the corresponding derived class.

.. code-block:: c++

    suffix = FileIO::getSuffix( fileName );
    std::unique_ptr<FileIO> file( FileIO::create( suffix );
    file->open( .. );
    ...

Distributed I/O Mode
====================

By default, a file is opened in the single mode. 

.. code-block:: c++

    if ( file->getDistributedIOMode() == DistributedIOMode::SINGLE )
    { 
        ...
    }

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

Extension for Other I/O Formats
-------------------------------

Adding support for any new file format is rather straightforward by writing
a class that derives from FileIO. 

.. code-block:: c++

    class MyIO : public FileIO,
                        FileIO::Register<MyIO> 
    {
        static std::string createValue();   // registration value for factory

        static FileIO* create();            // create routine called for create( createValue() )

        void writeStorage( const _MatrixStorage& storage );

        void readStorage( const _MatrixStorage& storage );

        void readArray( hmemo::HArray<ValueType>& array );
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


