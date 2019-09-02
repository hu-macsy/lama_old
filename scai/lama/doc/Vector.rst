.. _lama_Vector:

Vector
======

The class ``Vector`` is a generic mathematical vector and any object of this class stands for a one-dimensional array with values
of a certain type. In contrary to the ``HArray`` or ``LArray``, a ``Vector`` can be distributed among processors by having
a specific ``Distribution`` from :ref:`scaidmemo:main-page_dmemo`.
For a distributed vector, each processor stores only the values owned by it. The ``Vector`` class comes in two flavors:
a ``DenseVector`` stores all its values and a ``SparseVector`` stores only the non-zero values. In the latter case, an additional
array containing the indexes of the non-zero values is required, but the array containing the values might be significantly 
smaller. The local values are internally stored in a ``HArray`` out of :ref:`scaihmemo:main-page_hmemo` 
so a ``Vector`` can be used transparently on every device. 

Furthermore, each vector object has a context, that can be set explicitly. This context decides on which device operations 
with this vector are executed. In case of operations with multiple vectors it is usually the target vector that decides, where
an operation is executed.

Constructors
------------

For instantiating a vector variable you need to call the constructor either of the templated class ``DenseVector``, 
a specific representation of a vector holding all vector entries, or of the templated class ``SparseVector``, a
specific representation of a vector holding only non-zero entries.

.. code-block:: c++

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();
    const IndexType n = 1024;
    DistributionPtr blockDist( new BlockDistribution( n, dmemo::Communicator::getCommunicatorPtr() ) );
    ValueType initVal = 1;

    DenseVector<ValueType> vector1;
    SparseVector<ValueType> vector2( ctx );
    SparseVector<ValueType> vector3( n [, ctx] );
    DenseVector<ValueType> vector4( blockDist [, ctx] );
    SparseVector<ValueType> vector5( n, initVal [, ctx] );
    SparseVector<ValueType> vector6( blockDist, initVal [, ctx] );
    SparseVector<ValueType> vector7( vectorX );
    SparseVector<ValueType> vector8( vectorX, blockDist );
    SparseVector<ValueType> vector9( "vector.mtx" );            // read vector from a file
    DenseVector<ValueType> vector10( vector1 + 2 * vector3 );   // any supported vector expression

Sparse and dense vectors have syntactically nearly the same constructors.
Only the following constructors are available only for a sparse vector:

.. code-block:: c++

    HArray<IndexType> nonZeroIndexes;
    HArray<ValueType> nonZeroValues;
    ...   // compute set local non-zero values independently
    // Note: nonZeroIndexes.size() == nonZeroValues.size() must be valid
    SparseVector<ValueType> sVector1( n, nonZeroIndexes, nonZeroValues, zeroValue [, ctx] );

    IndexType indexes_raw[] = { 1, 3, 7 };
    ValueType values_raw[] = { 2, 1, 6 };
    SparseVector<ValueType>( n, 3, indexes_raw, values_raw [, initVal ] );

All index positions that do not appear in the array nonZeroIndexes are assumed to be zero. But keep
in mind that the zero element of a sparse vector might be any value and not necessarily the value 0.


The use of the following constructors is not really recommended as it results in vector that have
undefined values.

.. code-block:: c++

    SparseVector<ValueType> vector3( n );
    DenseVector<ValueType> vector4( blockDist );

Furthermore, keep in mind that the size of a vector or a distributon is set or changed by most operations
on vectors and so this attribute does not hold for the lifetime of a vector at all.

Vector Initialization
---------------------

Usually a vector should be allocated and fully initialized by its constructor. There might
be some situations where this is not possible:

 * A corresponding constructor is not available, e.g. for setting random numbers or setting
   a range of values like lb, lb + inc, lb + 2 * inc, ....
 * The vector might have been created dynamically where these routines usually generate 
   only a zero vector, i.e. not allocated and not initalized at all.

.. code-block:: c++

    DenseVector<ValueType> v1;

    v1.setData( array );          // initialize a vector with data from a heterogeneous array
    v1.setRawData( n, rawData );  // initialize a vector with any 'raw' data ( size, pointer )
    v1.setRange( n, 3, 2 );       // initializes the vector with the values 3, 5, 7, 9, ...
    v1.setRandom( n, 10 );        // initialize the vector with n random numbers in the range 0..10
    v1.setSameValue( n, 5 );      // initialize the vector with n elements of the value 5
    v1.readFromFile( "v.mtx" );   // read the vector from a file, size can queried afterwards

Sometimes it might be useful to create a sparse vector, i.e. a vector where most entries have the same
value, and only some values are different.

.. code-block:: c++

   v.setSparseData( n, zeroValues,  nonZeroIndexes, nonZeroValues );

Also a sparse random vector might be created where here a fill specifies the number of entries.

.. code-block:: c++

    v1.setSparseRandom( n, 0, 0.1f, 10 );   // initialize the vector with a certain ratio of random values

Most initialization routines might be called with a distribution instead of a size n. The initializations
of the local parts will be done independently.

.. code-block:: c++

    DenseVector<ValueType> v1;
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    v1.setRange( dist, 3, 2 );            // initializes the vector with the values 3, 5, 7, 9, ...
    v1.setRandom( dist, 10 );             // initialize the vector with n random numbers in the range 0..10
    v1.setSameValue( dist, 5 );           // initialize the vector with n elements of the value 5
    v1.setLocalData( dist, localArray );  // each processor initializes its local part with its data

In the first place it seems to be strange to use the initializon routines always with a size
or a distribution argument even if the vector has already been allocated. There are usually 
counterparts of these routines that do not require this first argument.

.. code-block:: c++

    DenseVector<ValueType> v( dist );   // define a distributed vector

    for ( int iter = 0; iter < MAX_ITER; ++iter )
    {
        v = 0;   // instead of v.setSameValue( dist, 0 );
        ...
    }

Nevertheless the use of the set routines are recommended for the following reasons:

* The size or distributon argument makes your code more stable and will even work
  if the vector has not been allocated or initalized before.
* There will be never any reallocation of memory as long as the size or distribution does not change.

Please not that for safety it is always a good strategy to initialize vectors with their allocation.
So in the following example code 1 might be more reliable than code 2 as in code 2 the allocated
vectors has undefined values between constructor and the call of the fill routine. But it is less efficient
as it does a complete write of the full vector data during the initializaton with 0 that is not 
required at all.

.. code-block:: c++

    // Code 1                          
    DenseVector<ValueType> v( n, 0 );   
    ...
    v = scalarValue;

    // Code 2
    DenseVector<ValueType> v( n );   
    ...
    v = scalarValue;    

    // Code 3
    DenseVector<ValueType> v;
    ...
    v.setSameValue( n, scalarValue ); 

The code 3 has the same efficiency as code 2 but it is more safe. This is due to the fact
that a zero vector causes less problems than an undefined allocated vector.

Vectors should be reused wherever it is possible. In the following loop the vector
is allocated and deallocated in each iteration of the loop even if the value n
is always the same.

.. code-block:: c++


    for ( int iter = 0; iter < MAX_ITER; ++iter )
    {
        DenseVector<ValueType> v( n, myValue );
        ....
    }

This code reuses the vector data in each iteration of the loop. Reallocation is
only done, if the value n becomes larger than any value used before.

.. code-block:: c++

    DenseVector<ValueType> v;

    for ( int iter = 0; iter < MAX_ITER; ++iter )
    {
        v.setSameValue( n, myValue );
        ....
    }

This is also one reason why you will never find any routine in LAMA that returns
a vector or a matrix. All supported vector operations in LAMA will never return a
new created vector. In the following example the implementation of the operator+ does not return
a vector but a syntactical construct that is resolved in the assignment and ends up in
an element-wise addition in-place in the exisiting vector v1. 

.. code-block:: c++

     DenseVector<ValueType> v1, v2, v3;
     ...
     v2.setXXX( ... )
     v3.setYYY( ... );
     v1 = v2 + v3;

Methods
-------

The class ``Vector`` is an abstract class that can be used for generic algorithm formulation. 
Beside some exceptions, all methods and vector expressions are supported for all kind of vectors,
either sparse or dense.

.. code-block:: c++

    Vector& v1 = denseVector;  
    Vector& v2 = sparseVector;

    v1.setContextPtr( Context::getContextPtr( Context::Host ) );
    v2.setContextPtr( Context::getContextPtr( Context::CUDA ) );

    const IndexType n = 100;
    v1.allocate( n );
    v2.allocate( DistributionPtr( new BlockDistribution( n, comm ) ) );
   
    v1 = ValueType( 2 );
    v2 = ValueType( 1 );

    v1.setDenseValues( denseValues );
    v1.setSparseValues( sparseIndexes, sparseValues, zeroValue );

    v2.readFromFile( "vector.mtx" );

For creating a new vector you need two major things:

 * the size of the vector (number of elements)
 * the value(s) of the vector

For distributed vectors the size can be substituted by a ``Distribution`` (holding the size and distribution strategy). 
For defining a Distribution, please refer to :ref:`this <scaidmemo:main-page_dmemo>` page.

The values can be passed by raw data pointer. Passing one value, will initilize the whole vector with this one value. 
Alternatively you can read the whole vector (size and data) from file, by specifing the filename. 
For a detailed description of the supported file formats, please refer to :ref:`lama_IO`.

Optionally you can specify a (initial) ``Context`` for the Vector, to define on which context the (initial) data is valid. 
For detailed explanation of the Context class, please refer to :ref:`this <scaihmemo:main-page_hmemo>` page. 

In the following you see all possible constructor calls:

.. code-block:: c++

  // for later use:
  const int size = 4;
  dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr( CommunicatorType::MPI ) );
  dmemo::DistributionPtr dist( dmemo::Distribution::getDistributionPtr( "BLOCK", comm, size, 1.0 ) );
  common::ContextPtr cudaContextPtr = common::Context::getContextPtr( common::context::CUDA );

  // empty (not initialized) float vector (with context, distribution, or both)
  DenseVector<float> empty();
  DenseVector<float> emptyDist( dist );
  DenseVector<float> emptyCUDA( cudaContextPtr );
  DenseVector<float> emptyDistCUDA( dist, cudaContextPtr );

  // creating a simple double Vector of size 4 with all elements having the value 1.0
  // optional third parameter: cudaContextPtr (hmemo::ContextPtr)
  DenseVector<double> x ( size, 1.0 );
  DenseVector<double> x2( dist, 1.0 );

  // reading from file (only on local vectors, can be redistributed afterwards)
  DenseVector<double> z( "z_vector.mtx" );

  // copy constructor (also works with general Vector 'z')
  DenseVector<double> zCopy   ( z );
  DenseVector<double> zRedist ( z, dist ); // z with a new Distribution

You also can create a pointer of a general Vector by calling the vector factory with a ``VectorCreateKeyType`` containing the vector type and the value type. The pointer can be saved as you need it as ``Vector*``, ``shared_ptr<Vector>``, ``unique_ptr<Vector>``. In LAMA we often make use of shared_ptr, so there is typedef to ``VectorPtr`` for that.

.. code-block:: c++

  // creating a DenseVector of value type double from the factory
  VectorCreateKey v_key( Vector::DENSE, common::getScalarType<double>() );
  VectorPtr vec_ptr = VectorPtr( Vector::create ( v_key ) );

For creating another Vector of the same type as your origin, you can receive the ``VectorCreateKeyType`` from it by calling ``getCreateValue()`` or ``getValueType`` for just getting the ValueType.

.. code-block:: c++

  VectorPtr z_clone1 = VectorPtr( Vector::create( z.getCreateValue() ) );              // or
  VectorPtr z_clone2 = VectorPtr( Vector::create( Vector::DENSE, z.getValueType() ) );

Vector Operations
------------------

Operations for sparse and dense vectors are the same as for LArrays.

.. code-block:: c++

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10;

    DenseVector<double> x( n, 1.0, ctx );
    DenseVector<double> y( n, 2.0, ctx );

    x[0] = 0.5;
    y[1] = x[0] * 1.0 - 0.5 * y[0];

    x += 1.0;
    y -= 1.3;
    y *= 1.5;
    x /= 0.7;

    x += y;
    y -= x;
    x /= y;
    x *= y;

    y += x *= 2;

    // unary operations

    x.invert();      // x[i] = 1.0 / x[i]
    y.conj();        // y[i] = conj( y[i] )
    x.log();
    y.floor();
    x.ceil();
    x.sqrt();
    x.sin();
    x.cos();
    x.tan();
    x.atan();
    x.powBase( 2.0 );  // x[i] = 2.0 ** x[i] 
    y.powExp( 2.0 );   // x[i] = x[i] ** 2.0
    x.powBase( y );    // x[i] = y[i] ** x[i]
    y.powExp( x );     // y[i] = y[i] ** x[i]

    Scalar s;

    s = x.sum();
    s = x.min();
    s = x.max();

    s = x.l1Norm();
    s = x.l2Norm();
    s = y.maxNorm();
   
    s = x.dotProduct( y );
    s = x.maxDiffNorm( y );

DenseVector or SparseVector
---------------------------

The following differences between a dense and a sparse vector should be kept in mind:

* There is no method to set individually a single element in sparse vector, while a dense vector has the method ``setValue``.
* gather and scatter operations are only supported for dense vectors
* sorting is only supported for dense vectors
* Many operations where vectors are involved require an explicit array with all (local) values. For a
  dense vector the method ``getLocalValues`` gives a reference to the corresponding heterogeneous array for free,
  for a sparse vector this array will be built temporarily by calling the method ``buildLocalValues``.

As a fallback, many methods use a dense array with all local values of a method. In these cases,
a sparse vector might perform slower than a dense vector. The following code shows the typical pattern
how to implement code that requires individual solutions, either if the vector is dense or sparse.

.. code-block:: c++

    const Vector& v = ...

    switch ( v.getVectorKind() )
    {
        case Vector::DENSE:
        {
            const _DenseVector& denseV = reinterpret_cast<const _DenseVector&>( v );
            ... denseV.getLocalValues()  ...  // only for dense vectors available
            break;
        }
        case Vector::SPARSE:
        {
            const _SparseVector& sparseV = reinterpret_cast<const _SparseVector&>( v );
            HArray<ValueType> v;
            sparseV.buildLocalValues( v );
            ...
            break;
        }
        default:
            COMMON_THROWEXCEPTION( "illegal vector kind: " << v.getVectorKind() )
    }

Here are some typical situtations where an application might benefit from a sparse vector:

- getRow or getColumn of a sparse matrix is faster if the result is stored in a sparse vector
- many binary operations with a dense and a sparse vector are faster, as shown in the following code

.. code-block:: c++

   Matrix& m;
   _SparseVector& sparseV = ...
   _DenseVector& denseV = ...

   m.getRow( sparseV, i );
   m.getColumn( sparseV, j );

   Scalar s = sparseV.dotProduct( denseV );
   Scalar s = denseV.dotProduct( sparseV );
   denseV += alpha * sparseV;
   denseV -= alpha * sparseV;

Binary operations with two sparse vectors (if not the same) require some overhead to determine the new pattern
for the non-zero elements.

.. code-block:: c++

   _SparseVector& sparseV1 = ...
   _SparseVector& sparseV2 = ...
   
   Scalar s = sparseV1.dotProduct( sparseV2 );
   sparseV1 += sparseV2;

Expressions
-----------

Having vectors and scalars (as ``Scalar`` or value) you can perform vector addition, substraction and scaling with a scalar in text-book syntax. We have implemented the expressions to a maximal length of the form:

.. code-block:: c++

    v_z = s_alpha * v_x + s_beta * v_y;

All specialization of this form (e.g. s_alpha = 1, s_beta = 0) are valid expressions:

.. code-block:: c++

    Scalar s( 2.0 );
    x = s * x;
    
    z = x + y;
    z = x * 2.0 + y;
    z = 2.0 * x + y;
    z = x + y * 1.0;
    
    z = y * 2.0;
    z = y / 2.0;
    
Also the combination with the assign operator is possible (internally handled as v_z = s_alpha * v_x + s_beta * v_z):

.. code-block:: c++

    z += x;
    z += 2.0 * x;
    z += x * 2.0;

    z -= x;
    z -= 2.0 * x;
    z -= x * 2.0;
    z *= 3.0;
    z /= 1.5;

For initializing a Vector, you can assign one value to the whole vector by the assignment operator ('='). The size of the vector is kept.

.. code-block:: c++

    x = 1.0;
    y = 2.0;

FFT
---

The following example shows how to call the Fast Fourier Transform for a vector (in-place):

.. code-block:: c++

   #include<scai/lama/fft.hpp>

   auto x = read<DenseVector<ComplexDouble>>( "input.mtx" );

   // Note: size of x must be a power of 2

   fft( x );    // apply fast fourier transform
   ifft( x );   // apply inverse fast fourier transform

   x.writeToFile( "output.mtx" );

Here are some remarks about calling fft or ifft for a vector:

 * The size of the vector must be a power of 2
 * The distribution does not change but it might be redistributed intermeadiately
 * The value type of the vector must be a complex type.


Utility Functions
-----------------

Additionally you have some utility functions that can be called on a vector: (for getting the size or distribution of the vector, e.g. after reading it from file, for swapping with another vector or creating a copy.

.. code-block:: c++

    IndexType length = x.size(); // getting the global size of a vector
    DistributionPtr d = x.getDistributionPtr(); 

    x.swap( y ); // swapping the size and values of the vectors

    Vector* zCopy = z.copy(); // calls the copy constructor

For accessing single values of a vector you can use ``getValue`` or ``()`` with the global index ``i``. But you must have in mind, that it may be inefficient if the vector is distributed and/or not on the Host Context, because of communication between nodes or CPU and GPU:

.. code-block:: c++

    s = z.getValue( index );
    s = z( index );

File I/O
--------

Except from a constructor with a passed string, you can use ``readFromFile`` and ``writeToFile``. The generally excepted format in LAMA for vector and matrices is defined :doc:`here<FileIO>`.

.. code-block:: c++

    x.readFromFile( "vector.mtx" );
    // writing a vector to file in matrix market format in double precision
    y.writeToFile( "result.mtx", File::MATRIX_MARKET, File::DOUBLE );

Math Functions
--------------

The dot product of two vectors is expressed as function ``dotProduct``:

.. code-block:: c++

    s = x.dotProduct( y );

Also the rudimental math functions 'max', 'min', are prepared on a ``Vector``, returning the global maximum/minimum of all entries.

.. code-block:: c++ 

   Scalar maximum = x.max();
   Scalar minimum = y.min();

You can get the L1-, L2-, Maximum-norm of an ``Vector`` by:
   
.. code-block:: c++ 
   
    s = x.l1Norm();
    s = x.l2Norm();
    s = x.maxNorm();

Output operator
---------------

Also the output operator for a ``Vector`` is implemented, giving you informations about its size, ``Distribution`` and ``Context``.

.. code-block:: c++ 
  
    std::cout << "my vector x looks like: " << x << std::endl;

The output will look like the following, telling you x is a DenseVector of type double with global and local size of four (therefore having a NoDistribution of size four that is located on the Host (CPU with 4 enabled OpenMP threads) ).

.. code-block:: bash

  my vector x looks like: DenseVector<double>( size = 4, local = 4, dist = NoDistribution( size = 4 ), loc  = HostContext( #Threads = 4 ) )
