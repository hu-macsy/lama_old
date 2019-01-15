.. _CollectiveFile:

CollectiveFile
==============

A collective file is a common file for all processors of a communicator from which it has been created.

.. code-block:: c++

    CommunicatorPtr commWorld = Communicator::getCommunicatorPtr();
    std::unique_ptr<CollectiveFile> file = commWorld->collectiveFile();
    // more convenient: auto file = commWorld->collectiveFile();

Each processor of the communicator can read data from or write data into this common file simultaneously.
An efficient implementation, e.g. using the MPI I/O implementation for MPI communicators,
allows for parallel I/O that can give better performance than serial or sequential I/O.

Methods of a collective file must be called simultaneously by all processors of the
communicator from which it has been created. While the method writeSingle writes 
replicated data (same values on all processors) once into a file, the method
writeAll writes distributed data (individual values on all processors), i.e.
the values of each processor in the file.

.. code-block:: c++

    file->open( "Test.data", "w" );
    HArray<ValueType> repArray = ...;   // same values on each processor
    HArray<ValueType> myArray = ...;   // individual values of each processor
    file->writeSingle( repArray );
    file->writeAll( myArray );
    file->close();

The following figure shows how the final collective file looks like.
In a similiar way data can be read from the file. It should be noted
that the number of processors that read the file can be different from
the number of processors that have written the file.

.. figure:: _images/collective_io.*
    :width: 500px
    :align: center
    :alt: CollectiveIO

    Collective File: write and read operations

Here is the corresponding code for reading a collective file.
The user has to take care that the number of read items (i.e. sizes
of the arrays) matches the number of entries in the file.

.. code-block:: c++

    file->open( "Test.data", "r" );
    HArray<ValueType> repArray;
    HArray<ValueType> myArray;
    file->readSingle( repArray, repSize );
    file->readAll( myArray, mySize );
    file->close();

Here are some remarks regarding the current implementation in LAMA:

 * Collective files are always considered to be binary files. 
 * read and write operations are only supported for those data types for which
   communication is supported.
 * writeAll and readAll write and read always contiguous sections for the 
   individual processors.
 * writeAll and readAll might exploit the possibilities of an underlying parallel
   file system.

In contrary to this collective or concurrent I/O, LAMA supports also independent
I/O where each processor writes its data into or reads it from an individual
file. This approach achieves usually better performance but has the big disadvantage
that the individual files can only be used by applications with the same number of
processors unless there is some tricky stuff to deal with other number of processors.
Independent I/O does not require any special features as it is just usual I/O 
individually for each processor and therefore out of the scope of the dmemo module.

.. figure:: _images/independent_io.*
    :width: 500px
    :align: center
    :alt: CollectiveIO

    Indepdendent I/O:
