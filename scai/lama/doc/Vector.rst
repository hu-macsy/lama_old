.. _lama_Vector:

Vector
======

The class ``Vector`` is a generic mathematical vector. The values are internally stored in a ``HArray`` out of :ref:`scaihmemo:main-page_hmemo` so a ``Vector`` can transparently used on every device. 

Constructors
------------

The class ``Vector`` is an abstract class that can be used for generic algorithm formulation.
For instantiating a vector variable you need to call the constructor of the templated class ``DenseVector``, a specific representation of a vector holding the full vector entries.

For creating a new DenseVector you need two major things:
 * the size of the vector (number of elements)
 * the value(s) of the vector

For distributed vectors the size can be substituted by a ``Distribution`` (holding the size and distribution strategy). For defining a Distribution, please refer to :ref:`this <scaidmemo:main-page_dmemo>` page.

The values can be passed by raw data pointer. Passing one value, will initilize the whole vector with this one value. 
Alternatively you can read the whole vector (size and data) from file, by specifing the filename. For a detailed description of the supported file formats, please refer to :doc:`IO`.

Optionally you can specify a (initial) ``Context`` for the Vector, to define on which context the (initial) data is valid. For detailed explanation of the Context class, please refer to :ref:`this <scaihmemo:main-page_hmemo>` page. 

In the following you see all possible constructor calls:

.. code-block:: c++

  // for later use:
  const int size = 4;
  dmemo::CommunicatorPtr comm( dmemo::Communicator::getCommunicatorPtr( "MPI" ) );
  dmemo::DistributionPtr dist( dmemo::Distribution::getDistribution( "BLOCK", comm, size, 1.0 ) );
  common::ContextPtr cudaContextPtr = common::Context::getContextPtr( common::context::CUDA );

  // empty (not initialized) float vector (with context, distribution, or both)
  DenseVector<float> empty();
  DenseVector<float> emptyDist( dist );
  DenseVector<float> emptyCUDA( cudaContextPtr );
  DenseVector<float> emptyDistCUDA( dist, cudaContextPtr );

  // creating a simple double Vector of size 4 with all elements having the value 1.0f
  // optional third parameter: cudaContextPtr (hmemo::ContextPtr)
  DenseVector<double> x ( size, 1.0f );
  DenseVector<double> x2( dist, 1.0f );

  // creating a local (not distributed) vector from raw double pointer
  const double inputData[] = { 1.0, 2.0, 3.0, 4.0 };
  scai::lamaDenseVector<double> y( size, inputData ); // optional third parameter: cudaContextPtr (hmemo::ContextPtr)

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

Except from a constructor with a passed string, you can use ``readFromFile`` and ``writeToFile``. The generally excepted format in LAMA for vector and matrices is defined :doc:`here<IO>`.

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
