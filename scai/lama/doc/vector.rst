.. _vector:

Vector
======

The class **Vector** is a generic mathematical vector. The values are internally stored in a HArray out of :ref:`hmemo <scaihmeo:main-page_hmemo>` so a Vector can transparently used on every device. 

Constructors
------------

The general class Vector is an abstract class that can be used for generic algorithm formulation.
For instantiating a vector variable you need to call the constructor of **DenseVector**, a specific representation of a vector holding the full vector entries.

For creating a new DenseVector you need two major things:
 * the size of the vector (number of elements)
 * the value(s) of the vector

For distributed vectors the size can be substituted by a **Distribution** (holding the size and distribution strategy). For defining a Distribution, please refer to :ref:`this <scaidmemo:main-page_dmemo>` page.

The values can be passed by raw data pointer. Passing one value, will initilize the whole vector with this one value. 
Alternatively you can read the whole vector (size and data) from file, by specifing the filename.

Optionally you can specify a (initial) **Context** for the Vector, to define on which context the (initial) data is valid. For detailed explanation of the Context class, please refer to :ref:`this <scaihmemo:main-page_hmemo>` page. 

In the following you see all possible constructor calls:

.. code-block:: c++

  // creating a simple float Vector of size 8 with all elements having the value 1.0f
  DenseVector<float> x( 8, 1.0f );

  // creating a local (not distributed) vector from raw double pointer
  const int size = 4;
  const double inputData[] = { 1.0, 2.0, 3.0, 4.0 };
  DenseVector<double> y( size, inputData ); // optional third parameter: hmemo::ContextPtr

  //alternatively: the distributed versions
  DistributionPtr dist ( ... );
  DenseVector<float> x2( dist, 1.0f );
  DenseVector<double> y2( dist, inputData ); // optional third parameter:: hmemo::ContextPtr

  // reading from file
  DenseVector<float> z( "in.txt" );

  // 


Having vectors and scalars (as Scalar or value) you can perform vector addition, substraction and scaling with a scalar
in text-book syntax:

.. code-block:: c++

    Scalar s( 2.0 );
    x = s * x;
    
    z = x + y;
    z = x * 2.0 + y;
    z = 2.0 * x + y;
    z = x + 1.0 * y;
    z = x + y * 1.0;
    
    z = y * 2.0;
    z = y / 2.0;
    
and also combined with the assign operator:

.. code-block:: c++

    z += x;
    z += 2.0 * x;
    z += x * 2.0;

    z -= x;
    z -= 2.0 * x;
    z -= x * 2.0;
    z *= 3.0;
    z /= 1.5;

You can assign one value to the whole vector also by '=' (size is kept):

.. code-block:: c++

    x = 1.0;
    y = 2.0;
    
The dot product of two vectors is expressed as function (so there is no misunderstanding):

.. code-block:: c++

   s = x.dotProduct( y );

You can get the L1-, L2-, maximum-norm of an vector by:
   
.. code-block:: c++ 
   
   s = x.l1Norm();
   s = x.l2Norm();
   s = x.maxNorm();

Other useful functions on a vector are:

.. code-block:: c++ 

   IndexType length = x.size(); // getting the size of a vector
   
   Scalar maximum = x.max(); // getting the maximum value of all entries
   Scalar minimum = x.min(); // getting the minimum value of all entries
   
   // writing a vector to file as formatted output in double precision
   z.writeToFile( "vector.txt", File::FORMATTED, File::DOUBLE);