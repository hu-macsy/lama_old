.. _lama_Scalar:

Scalar
======

The class ``Scalar`` represents a scalar value, storing every value type. It can be a value in an expression, or the result of an expression. Compared to ``Vector`` and ``Matrix`` a ``Scalar`` can not be distributed amoung nodes nor is located on a specific location (it is always created on the host and will be passed to other locations).

Constructors
------------

For initializing a ``Scalar`` call the constructor.
You can create a ``Scalar`` by delivering a value to the constructor. The value type delivered is preserved for the scalar.
An empty constructor prepares a scalar representing zero. Delivering another ``Scalar`` copies the value, preserving the value type.

.. code-block:: c++

   Scalar a( 1 );                         // a scalar of IndexType (int)
   Scalar b( 2.5f );                      // a scalar of type float
   Scalar c( 3.14 );                      // a scalar of type double
   Scalar d( ComplexDouble( 2.0, 1.5 ) ); // a scalar of type ComplexDouble
   Scalar e();                            // a scalar representing zero
   Scalar f( c );                         // a new scalar with the value and value type of c

Expressions
-----------

For calculating with a ``Scalar`` all common binary operators '+', '-', '\*', '/' (also as '+=', '-=', '\*=','/=') and the unary operator '-' are supported. A mix of different template types is allowed.

.. code-block:: c++

   // binary operators
   c = a + b;
   c = a - b;
   c = a * b;
   c = a / b;
   
   a += b;
   a -= b;
   a *= b;
   a /= b;
   
   // unary operator '-'
   c = -c;

Math functions
--------------

Also the rudimental math functions 'max', 'min', 'abs', 'sqrt' are prepared and can be used as known from the std.

.. code-block:: c++

   c = max( a, b );
   c = min( a, b );
   c = abs( a );
   c = sqrt( a );


Comparation operators
---------------------

All known relational operators '==', '!=', '<', '>', '<=', '>=' are overloaded and can be used as the following:

.. code-block:: c++

   bool boolean = ( a == b );
   boolean = ( a != b );
   boolean = ( a < b );
   boolean = ( a > b );
   boolean = ( a <= b );
   boolean = ( a >= b );

Output operator
---------------

For printing a ``Scalar`` value you can use the output operator ('<<'), resulting in an output of the form ``Scalar(<your_value>)``.

.. code-block:: c++

   std::cout << "my Scalar is: " << a << std::endl;

The output will be the following:

.. code-block:: c++

   my Scalar is: Scalar(1)

Alternatively you can print the plain value on your own, e.g. in printf(). You receive the value in the desired value type by calling 'getValue<ValueType>()'

.. code-block:: c++
   
   // getValue
   printf("my Scalar is: %d\n", a.getValue<int>() );
   