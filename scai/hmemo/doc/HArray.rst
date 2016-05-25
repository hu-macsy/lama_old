.. _HArray:

Heterogeneous Array
===================

HArray is a template class where the template argument specifies the type of the elements of the array.

.. code-block:: c++

    HArray<double> A( N, 1.0 );

Compared to other container classes in C++ like the vector class, the following differences
should be observed:

 * Default constructor of ValueType is not exploited, explicit initialization is necessary.
 * Copy constructor of ValueType is not exploited, memory transfer is just done bytewise.

Usually the value type is a primitive type or other arithmetic types like complex<float> or complex<double>.

HArray
------

HArray is container where data might be allocated at different memory locations. 

One incarnation is an object of the class HData and contains the following member variables:

=========    ==============================================
Property     Value
=========    ==============================================
mMemory      memory location 
mPtr         pointer to the allocated data
mCapacity    number of allocated bytes
mValid       true if this incarnation contains valid values
=========    ==============================================

Member variables of HArray:

==========   ==============================================
Property     Value
==========   ==============================================
mSize        size of the array
mValueSize   number of bytes for one element
constFlag    true for read-only
mData        list of current incarnations
mAccess      Read/Write access 
==========   ==============================================

Snapshot of a HArray( size = 1024, valueSize = 8 )

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
CUDA-0   8192         false
CUDA-1   8192         false
=======  ==========   ========

After initialization, at least one data entry is always valid.
The capacity specifies the allocated size on the memory, must be sufficient in case of valid data.

Here are some examples of initialization:

.. code-block:: c++

  HArray<double> A( 1024 );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
x        x            x
=======  ==========   ========

.. code-block:: c++

  HArray<double> A( host );

A( size = 0, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     0            false
=======  ==========   ========

.. code-block:: c++

  HArray<double> A( 1024, host );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         false
=======  ==========   ========

.. code-block:: c++

  HArray<double> A( 1024, host, 1.0 );

A( size = 1024, valueSize = 8 ):

=======  ==========   ========
memory   capacity     valid
=======  ==========   ========
Host     8192         true
=======  ==========   ========


Alias Problem
-------------

Many mathematical routines might use lhs array also on rhs.

.. code-block:: c++

  template<typename T>
  void add ( HArray<T>& res, const HArray<T>& a, const Harray<T>& b )
  {
      SCAI_ASSERT_LE( a.size(), b.size(), "size mismatch" )
  
      IndexType n = a.size();

      ContextPtr hostCtx = Context::getHostPtr();

      WriteOnlyAccess<T> write( res, hostCtx, n );
      ReadAccess<T>readA( a, hostCtx );
      ReadAccess<T>readB( b, hostCtx );
   
      double* resPtr = write.get();
      const double* aPtr = readA.get();
      const double* bPtr = readB.get();
  
      add[context]( resPtr, aPtr, bPtr, n );
  }
 
.. code-block:: c++

  add( a, b, c ); // this is okay
  add( a, a, b ); // here we have an ALIAS

Solution 1:

Check for alias and use appropriate accesses:

.. code-block:: c++

  if ( &res == &a )
  {
      hmemo::WriteAccess<T> write( res, hostCtx, n );
      hmemo::ReadAccess<T>readB( b, hostCtx );
      add[context]( write.get(), write.get(), readB.get(), n );
  }
  else if ( &res == &b )
  {
      hmemo::WriteAccess<T> write( res, hostCtx, n );
      hmemo::ReadAccess<T>readA( a, hostCtx );
      add[context]( write.get(), readA.get(), write.get(), n );
  }
  else
  {
      hmemo::WriteOnlyAccess<T> write( res, hostCtx, n );
      hmemo::ReadAccess<T>readA( a, hostCtx );
      hmemo::ReadAccess<T>readB( b, hostCtx );
      add[context]( write.get(), readA.get(), readB.get(), n );
  }
 

Solution 2:

Allow write and read access at same context.

.. code-block:: c++

  hmemo::WriteOnlyAccess<T> write( res, ctx, n );
  hmemo::ReadAccess<T>readA( a, ctx);
  hmemo::ReadAccess<T>readB( b, ctx);

If res == a, then ReadAccess after WriteAccess is not allowed as data might be resized.

.. code-block:: c++

  hmemo::ReadAccess<T>readA( a, ctx);
  hmemo::ReadAccess<T>readB( b, ctx);
  hmemo::WriteOnlyAccess<T> write( res, ctx, n );

This is fine, resize on res is not needed.

Prefetch
--------

Each HArray has a prefetch method in order to get a valid incarnation of the array
at a certain context. If a memory transfer is required, this memory transfer is started
asynchronously. Synchronization is done implicitly with the next access to the array.

.. code-block:: c++

  ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA );
  ContextPtr hostContext = Context::getContextPtr( common::context::Host );
  HArray<double> A;
  ...
  {
      WriteAccess<double> wA( A, hostContext );
      ...
  } // valid data only on Host

  A.prefetch( cudaContext ); // starts async transfer Host->GPU

  workload( dummy, NWORK );  // overlaps with memory transfer

  {
      WriteAccess<double> rA( A, cudaContext ); // waits until transfer is complete
      ...
  }

An asynchronous memory transfer to a CUDA device is done via a CUDA stream if the 
incarnation of the array on the Host is in the CUDA Host memory. Otherwise a
separate thread is started that takes care of the memory transfer.

Using Pinned Memory
-------------------

By default, an incarnation of an Harray on the host is allocated in the Host memory.

.. code-block:: c++

  ContextPtr hostContext = Context::getContextPtr( common::context::Host );
  HArray<double> A;
  {
      // will use HostMemory
      WriteOnlyAccess<double> wA( A, hostContext, N );
      ...
  }

If the data is later needed on the GPU, no fast memory transfer is possible as 
the data is not in the pinned memory. The allocation of host memory in the pinned 
memory can be forced as follows:

.. code-block:: c++

  ContextPtr hostContext = Context::getContextPtr( common::context::Host );
  ContextPtr cudaContext = Context::getContextPtr( common::context::CUDA );

  HArray<double> A( cudaContext );
  {
      // will use CUDAHostMemory
      WriteOnlyAccess<double> wA( A, hostContext, N );
      ...
  }
  {
      ReadAccess<double> ra( A, cudaContext );   // fast memory transfer
      ...
  }

The use of a context pointer in the constructor works like a first touch on the 
corresponding context.

.. code-block:: c++

  HArray<double> A( cudaContext );

  HArray<double> A;
  {
      WriteAccess<double> dummyW( A, cudaContext );
  }

Some other strategies have been dropped for these reasons:

 * Using pinned memory as default memory on the Host is not recommended as 
   allocation in pinned memory is 10 up to 100 times slower.

 * Pinning already allocated unpinned memory might be possible e.g. when data transfer
   to the GPU is required. This does not fit well in the design concept that handles
   the two memory classes separately. Furthermore, using this as a general strategy
   is not always recommended as the overhead does not pay off with one single transfer.

HArrayRef
---------

Each incarnation of an HArray is allocated in a corresponding memory object where the
memory management is handled by a corresponding manager.
Therefore data from any input array must be copied explicitly in the container.

.. code-block:: c++

  double* data = new double[N];
  read_data( data, N );
  ...
  HArray<double> vector;
  {
      WriteOnlyAccess<double> write( vector, host, N );
      for ( IndexType i = 0; i < N; ++i ) write[i] = data[i];
  }

  ReadAccess<double> write( vector, gpu); // now work on it on GPU

The class HArrayRef is provided to deal with such a situaton.

.. code-block:: c++

  double* data = new double[N];
  read_data( data, N );
  HArrayRef<double> vector( data, N )
  WriteAccess<double> write( vector, gpu); // now work on it on GPU

The memory at the pointer data will be used as incarnation on the Host memory 
As the data is not copied, it is not possible to resize the HArray ``vector``. 
The destructor of the array takes care that the specified memory will contain
a valid copy of the data.

Non-Zero Copy
-------------

 * CUDA devices can also operate on CUDA Host memory
 * Avoids data transfer to and from the device
 * But: operations are much faster if CUDA device memory is used

.. code-block:: c++

  {
      // Init on host
      WriteOnlyAccess<double> write( array, hostContext, N );
      double* v = write.get();
      for ( IndexType i = 0; i < N; ++i )
      {
          v[i] = 0.0;
      }
  }
  {    // work on GPU
       WriteAccess<double> write( array, cudaContext );
       work1[cudaContext]( write.get(), N );
  }

  {    // work on Host
       WriteAccess<double> write( array, hostContext );
       work2[hostContext]( write.get(), N );
  }

If array is allocated in CUDAHostMemory, no data transfer is needed.
