.. _lama_DenseVector:

DenseVector
===========

The class ``DenseVector`` is a container class and any object of this class stands for a one-dimensional array with values
of a certain type. In contrary to the ``HArray``, a ``DenseVector`` is distributed among processors by having
a specific ``Distribution`` (see :ref:`scaidmemo:main-page_dmemo`).
For a distributed vector, each processor stores only the values owned by it. 
These owned local values are internally stored in a ``HArray`` (see :ref:`scaihmemo:main-page_hmemo`)
so a ``DenseVector`` can be used transparently on every device. 

.. figure:: _images/dense_vector.*
    :width: 600px
    :align: center
  
    Distributed dense vector with the owned local values.

Constructors
------------

The most essential constructor for the template class ``DenseVector``
takes a distribution and the array with the local values to build the dense
vector.

.. code-block:: c++

    const IndexType size = 980;
    DistributionPtr blockDist =  blockDistribution( size );
    HArray<ValueType> localValues = ...;   //   must have blockDist->getLocalSize() entries

    DenseVector<ValueType> v( blockDist, localValues[, ctx] );

The context argument of a vector decides where operations on this vector will be executed. It
is optional, by default the context is the currently actual context (as specified by the environment
variable ``SCAI_CONTEXT``).

Instead of the array with the local values a constant single value can be used to construct
a dense vector.

.. code-block:: c++

    const ValueType initVal = 1;
    DenseVector<ValueType> v( blockDist, initVal[, ctx] );

The default constructor creates an empty vector. It is provided to allow this class to 
be used in C++ container classes.

.. code-block:: c++

    DenseVector<ValueType> v;
    std::vector<DenseVector<ValueType>> manyV;

The use of the following constructor is not really recommended as it results in a vector that has
undefined values. Nevertheless it might be useful in applications where allocation is separated from 
initialization of data.

.. code-block:: c++

    DenseVector<ValueType> v( blockDist );   // values of remain undefined
    v = 1;                                   // now v is initialized


Intialization of DenseVector
----------------------------

A dense vector can be initialized by different methods. Here are some methods that
allocate and initialize the vector:

.. code-block:: c++

    DenseVector<ValueType> v;
    v.setSameValue( dist, 5 );
    v.setLocalData( dist, localArray );
    v.readFromFile( file );

The following methods set or reset the values of an already allocated vector.

.. code-block:: c++

     DenseVector<ValueType> v( dist );
     v.fillLinearValues( start, inc );
     v.fillRandom( bound );
     auto lambda = []( IndexType i ) { return ValueType( i ); };
     v.fillByFunction( lambda );

Free Functions with DenseVector Result
--------------------------------------

For convenience, many methods for the initializaton of vectors are also available
as free functions that return a corresponding dense vector.

.. code-block:: c++

    const auto vector1 = denseVectorLinear<ValueType>( dist, firstVal, incVal );
    const auto vector2 = denseVectorZero<ValueType>( dist );
    const auto vector2 = denseVector<ValueType>( dist, initVal );
    const auto vector3 = denseVectorEval( matrix * x );
    const auto vector4 = denseVectorRead<ValueType>( "matrix.frm" )
    const auto vector5 = denseVectorRead<ValueType>( file )

.. code-block:: c++

    template<typename ValueType>                  template<typename ValueType>
    DenseVector<ValueType> freeFunction( ... )    void subroutine( DenseVector<ValueType>& v, ... )
    {                                             {
        DenseVector<ValueType> v( ... );               v = ...;
        ...                                            ...
        return v;                                 }
    }

Dense Vector Operations
-----------------------

DenseVector is a derived class from the generic class ``Vector``, so all methods and 
operations provided by this class are also available for the ``DenseVector`` class.
This includes especially all vector operations.

.. code-block:: c++

    DenseVector<ValueType> v1 = ...;
    DenseVector<ValueType> v2 = ...;
    DenseVector<ValueType> v3 = ...;
    const Matrix<ValueType>& m  = ...;
    ValueType alpha = ...;
    ValueType beta  = ...;
    v1 = alpha * m * v2 + beta * v3;
    v1 = alpha * v2 + beta * v3;
    v1 = pow( v2, v3 );
    v1 = v2 / v3;
    alpha = v1.l1Norm();
    beta = v2.dotProduct( v3 );
    ...

For convenience, the index operator [] can be used to access single elements, but
is not recommeded.

.. code-block:: c++

    IndexType size = atoi( argv[2] );
    auto dist = blockDistribution( size );
    DenseVector<double> v( dist, 0.0 );
    v[0] = v[size-1] = 1.0;
    ...
    double m = v[ size / 2 ];

For implementation of own functions it is recommended to use either write or 
read accesses to the local elements.

.. code-block:: c++

    const DenseVector<double> v = ...;
    DenseVector<double> w;
    auto dist = v.getDistributionPtr();
    w.allocate( dist );
    auto writeW = hostWriteAccess( w.getLocalValues() );
    auto readV = hostReadAccess( v.getLocalValues() );
    for ( IndexType i = 0; i < dist->getLocalSize(); ++i )
    {
        writeW[i] = f( readV[i], ... );
    }

FFT
---

The following example shows how to call the Fast Fourier Transform for a vector (in-place):

.. code-block:: c++

   #include<scai/lama/fft.hpp>

   auto x = denseVectorRead<ComplexDouble>( "input.mtx" );

   // Note: size of x must be a power of 2

   fft( x );    // apply fast fourier transform
   ifft( x );   // apply inverse fast fourier transform

   x.writeToFile( "output.mtx" );

Here are some remarks about calling fft or ifft for a vector:

 * The size of the vector must be a power of 2
 * The distribution does not change but it might be redistributed intermeadiately
 * The value type of the vector must be a complex type.

Replicated DenseVector
----------------------

You might also create a replicated dense vector, i.e. a full dense vector, that has an incarnation
on each processor. Therefore you can use the class ``NoDistribution`` for specifying the distribution.

.. code-block:: c++

    const IndexType size = 1024;
    DistributionPtr noDist = noDistribution( size );
    ValueType initVal = 1;
    HArray<ValueType> values = ...;    // must have 'size' entries;
    DenseVector<ValueType> v1( noDist, initVal );
    DenseVector<ValueType> v2( noDist, values ); 

Most methods and functions for dense vectors allow the use of the size argument directly instead of the 
distribution (pointer) argument. If the size can be deduced from other arguments,
it does not have to be used at all.

.. code-block:: c++

    const IndexType size = 1024;
    ValueType initVal = 1;
    HArray<ValueType> values = ...;    // must have 'size' entries;
    DenseVector<ValueType> v1( size, initVal );
    DenseVector<ValueType> v2( values );   // size is same as values.size()
    const auto v3 = denseVectorZero<ValueType>( size );


Good Practice Advices
---------------------

A vector is not a container class where elements can easily be added or removed. Therfore other C++
container classes should be used, and a LAMA vector should only be generated when its size is known.

Like for all container classes, a dense vector allocates memory when it is constructed.
This allocation is not cheap at all and therefore vectors should be reused wherever possible.
Especially reuse in loops should have rather high priority in order to achieve good performance.

.. code-block:: c++

    auto x = denseVector( blockDistribution( N ), ValueType( 0 ) );

    // bad practice                                 // good practice
    for ( IndexType i = 0; i < NITER; ++i )         DenseVector<ValueType> v;
    {                                               for ( IndexType i = 0; i < NITER; ++i )
        auto v = denseVectorEval( matrix * x );     {
        x = x + alpha * v;                              v = matrix * x;
        ...                                             x = x + alpha * v;
    }                                                   ....
                                                    }

Even if the size or distribution of a vector might change during its lifetime, it is recommended
to avoid these operations, i.e. a dense vector should be allocated once with its size as used
for the application.

