.. _distributed_IO:

Distributed I/O Support in LAMA
===============================

In contrary to a serial program where all data of vectors and/or matrices reside  
in the local memory of one processor, the I/O of distributed data involves some
kind of management to deal with it.

Single I/O
----------

In the single I/O mode, only the first or master processor is responsible for
reading data from and writing data into the file. Therefore data will be
redistributed corresponding to the ``SingleDistribution`` where the 
first processor owns all data.

.. figure:: _images/single_io.*
    :width: 600px
    :align: center
  
    Single I/O of Distributed Data.

Advantages:

 * This mode supports all I/O formats that are supported in serial programs.
 * All data is available in a single file.

Disadvantages:

 * Data must be redistributed (involves communication)
 * First processor has to allocate memory for all data
 * Poor performance.

Independent I/O
---------------

.. figure:: _images/independent_io.*
    :width: 600px
    :align: center
  
    Independent I/O of Distributed Data.

Advantages:

 * This mode supports all I/O formats that are supported in serial programs.
 * Full parallel performance
 * Might use local discs (e.g. fast SSD) of the nodes

Disadvantages:

 * Lots of files to manage
 * Number of processors must match the number of files (reading).
 * Might involve redistribution if data is not block-distributed

When a read or write operation is applied to a *distributed* vector or matrix, it
is just as if the local part is read or written. But as the distribution itself is
not written into the file, the vector or matrix must be block distributed.
This implies a redistribution of the data before a write operation. The read data
will have a general block distribution where the local size corresponds to the
size of the data in the local files.

Note: when writing data with coordinates (vector in sparse format) or COO storage
data, the local files still keep the local coordinates. It might be possible to use
global coordinates and to avoid the redistribution before the write operation.

Collective I/O
---------------

In the collective mode, all processors read from the same 
file or write into it. The module ``dmemo`` provides the ``CollectiveFile`` 
class that can be used to read/write distributed data structures
like vectors and matrices.
 
.. figure:: _images/collective_io.*
    :width: 600px
    :align: center
  
    Collective I/O of Distributed Data.

Advantages:

 * All data resides in one file.
 * Number of processors for read can be different from number of processors for write.
 * Parallel file systems allow efficient implementation of read/write operations
   for collective files.

Disadvantages:

 * Only possible for file formats where processors can read/write from individual offsets.
 * Might involve redistribution if data is not block-distributed

Explain: what happens if data is not block distributed.

Explain: what happens with unsorted COO data. So we assume always sorted COO data.

Selection of Distributed I/O-Mode
---------------------------------

.. code-block:: c++

    #include <scai/lama/io/LamaIO.hpp>

    using namespace scai;

    LamaIO file;

    file.open( "data.lfm", "r", DistributedIOMode::SERIAL );
    if ( comm->getRank() == 0)
    {
        file.read( localArray );
    }
    else
    {
        localArray.clear();
    }
    DenseVector<double> v( localArray, 
    CSRSparseMatrix<double> m = read<CSRSparseMatrix<double>>( file );
    file.close();

    file.open( "vector.lfm", "w", DistributedIOMode::COLLECTIVE );
    v.writeToFile( file );
    m.writeToFile( file );
    file.close();

Partitioned I/O
================

The FileIO classes are used to read and write local data, i.e. heterogeneous arrays and matrix storages.
The Partioned I/O supports read and write of distributed data, i.e. vectors, matrices, and also distributions.

Two modes are supported for distributed data structures:

- Single File Mode: all data is written into one single file and read from it. In this mode the data
  is communicated to the master process that takes care of the file access.

- Multiple File Mode: in this mode each processor writes its local data into a separate file.

The single file mode is the default mode for the read and write method of distributed data.
The multiple file mode is chosen if the pattern "%r" appears in the file name. This pattern 
will be replaced with the corresponding rank and size values of the communcator.

.. code-block:: c++

    Vector& v = ...                            Vector& v = ...
    Matrix& m = ...                            Matrix& m = ...
    ...
    v.writeToFile( "vector.mtx" )              v.readFromFile( "vector.mtx" )
    m.writeToFile( "matrix.mtx" )              m.readFromFile( "matrix.mtx" )
    ...
    v.writeToFile( "vector%r.mtx" )            v.readFromFile( "vector%r.mtx" )
    m.writeToFile( "matrix%r.mtx" )            m.readFromFile( "matrix%r.mtx" )

The single file mode has the big disadvantage that all data is written and read by a single processor
and therefore memory is allocated for the full matrix by this master process. This might cause serious problems if 
matrices become too large. 

The multiple file mode causes less memory overhead and might be much faster especially on hardware
that supports parallel I/O. Currently it has the disadvantage that the number of processors for reading
distributed data must be exactly the same as the number of processors that have written the data.

When writing distributed data into a single or multiple file, information about the distribution itself
is lost. For a single file this does not really matter as the data can be redistributed arbirtrarily after reading
it. When the data is written into multiple files, the information about the mapping is lost, i.e. it is no more known
how the submatrix of one file fits into the whole matrix. Only for block or general block distributions the original
mapping from local to global indexes can be determined. Otherwise, the mapping itself must also be written into
a file or, as it is a distributed data structure, into multiple files.

Single File IO
--------------

Idea: all data is gathered on the master process and written to a single file.
When reading a matrix or vector from a single file, the master process reads it.

There are some pitfalls.

* Writing a replicated vector or matrix: it must be called by all processors even if only one
  processor writes it. It contains an implicit synchronization that avoids that the same file
  is read before it is completely written.

.. code-block:: c++

    Vector& v = ...;
    v.redistribute( DistributionPtr( new NoDistribution( v.size() ) ) );
    v.writeToFile( "vector.mtx" )          // implicit synchronization 
    v.readFromFile( "vector.mtx" )

General rule: each I/O routine of a Vecotor or a Matrix must be called by all processors
of the current communicator.

Reading a vector or matrix from a single file comes up with a distribution where only
the master process has the data.

.. code-block:: c++

   CylcicDistribution dist( globalSize, globalSize, comm );    // first processor owns all elements

I/O of Distributions
--------------------

If a vector or a matrix is stored in a single or multiple files, the information about the mapping is lost 
if the distribution is not a (general) block distribution.

Therefore the mapping itself can be written to a single or mutliple files. 

A single file contains for each entry the owner.

.. code-block:: c++

   PartitionIO::write( distribution, "owners.mtx" )

   CommunicatorPtr comm = Communcator::CommuncatorPtr();    // current MPI comm
   DistributionPtr dist( PartitionIO::readDistribution( "owners.mtx", comm )

.. code-block:: c++

   0   0   1    1   0   0   2   2   3   3   2   2   3   3

A multiple file contains for each processor the owned global indexes.

.. code-block:: c++

   0 :   0  1   4   5
   1 :   2  3   6   7
   2 :   8  9  12  13
   3 :  10 11  14  15

.. code-block:: c++

   PartitionIO::write( distribution, "owners%r.mtx" )

   CommunicatorPtr comm = Communcator::CommuncatorPtr();    // current MPI comm
   DistributionPtr dist( PartitionIO::readDistribution( "owners%r.mtx", comm )

Partitioned I/O of Vectors
--------------------------

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

I/O of Matrices
---------------

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

I/O Methods for Vector and Matrix classes
-----------------------------------------

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

To write a ``Matrix`` or ``Vector`` to file call *writeToFile* with the name of the output file and the formatting. The default for just giving a name and no formatting is binary SAMG format in internal precision for the *values* und int for *ia* and *ja*.

.. code-block:: c++

   csrMatrix.writeToFile( "matrix_out.mtx", File::MATRIX_MARKET, File::FLOAT );
   vec.writeToFile( "vec_out.frv", File::SAMG_FORMAT, File::DOUBLE, /*binary*/ true ); // binary SAMG format
   
Possible file formats are ``File::SAMG_FORMAT`` and ``File::MATRIX_MARKET``.

Possible data types are ``common::scalar::INDEX_TYPE`` (int), ``common::scalar::FLOAT``, ``common::scalar::DOUBLE``, ``common::scalar::COMPLEX``(ComplexFloat), ``common::scalar::DOUBLE_COMPLEX``, ``common::scalar::LONG_DOUBLE_COMPLEX`` or ``common::scalar::INTERNAL`` for the internal representation of the data.

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

