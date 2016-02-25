Read-/WriteAccess
-----------------

Each use of the heterogeneous array requires an access

 * Accesses take care of valid data at the given context
 * i.e. take care of the corresponding memory transfer
 * Check for consistent use ( no read/write accesses at different locations)
 * Get operation provides pointer to the data
 * Be careful: never use data pointers outside the scope of accesses
 * Most checks are only possible at runtime

Example of ReadAccess:

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
=======  ==========   ========

.. code-block:: c++

  ReadAccess<double> readA( A, cuda )

By the read access, memory is allocated at CUDA device memory
and valid data is copied from Host to CUDA.

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
CUDA     8192         true
=======  ==========   ========

Example of WriteAccess:

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
=======  ==========   ========

.. code-block:: c++

  WriteAccess<double> writeA( A, cuda )

By the write access, memory is allocated at CUDA device memory
and valid data is copied from Host to CUDA. As data will be written
on the CUDA device, the Host incarnation becomes invalid.

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         false
CUDA     8192         true
=======  ==========   ========

Example of WriteOnlyAccess:

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
=======  ==========   ========

.. code-block:: c++

  WriteOnlyAccess<double> writeA( A, cuda )

By the write access, memory is allocated at CUDA device memory
but **no data** is copied from Host to CUDA. As data will be written
on the CUDA device, the Host incarnation becomes invalid.

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         false
CUDA     8192         true
=======  ==========   ========

Example of resize:

A( size = 1024, valueSize = 8 )

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
CUDA-1   8192         false
CUDA-2   16384        true
=======  ==========   ========

.. code-block:: c++

  A.resize( 2048 );

The resize operation only reallocates data for all valid incarnations that
have not sufficient capacity. So for the above example only the Host
incarnation is reallocated.

A( size = 2048, valueSize = 8 )

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     16384        true
CUDA-1   8192         false
CUDA-2   16384        true
=======  ==========   ========

Clear, Purge:

clear is the same as resize( 0 ) operation. It does not free any memory.
If the size is 0, the valid flag of an incarnation does not matter.

.. code-block:: c++

  A.clear();

A( size = 0, valueSize = 8 )

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     16384        true
CUDA-1   8192         false
CUDA-2   16384        true
=======  ==========   ========

Access Conflicts
----------------

The following situations result in an access conflict:

 * Two accesses on a different context if at least one is a write access;
   there is no guarantee for valid data
 * Two accesses by different threads if at least one is a write access
   There is not guarantee for valid data

Read and write access by same thread on same context is possible
Due to alias (e.g. X = 5 * X + 3 * Y ) this happens. But with a write and 
read access at the same time, the resize operation will throw an exception
if a reallocation becomes necessary.

Example 1:

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );
  WriteAccess<double> writeA( A, cuda );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         false
CUDA     8192         true
=======  ==========   ========

.. code-block:: c++

  ReadAccess<double> readA( A, host);  // fails, conflict to write on CUDA
  WriteAccess<double> writeA1( A, host);  // fails, conflict to write on CUDA
  ReadAccess<double> readA( A, cuda);  // fails, as might have been resized 

Example 2:

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );
  ReadAccess<double> readA( A, cuda );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
CUDA     8192         true
=======  ==========   ========

.. code-block:: c++

  ReadAccess<double> readA1( A, host);  // okay
  WriteAccess<double> writeA( A, host);  // fails
  WriteAccess<double> writeA( A, cuda);  // okay for same thread
  writeA.resize( 2048 );                 // fails
