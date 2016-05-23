Scalar
======

The class **Scalar** represents a scalar value, storing every value type. It can be a value in an expression, or the result of an expression.

Constructors
------------

For initializing a **Scalar** call the constructor.
You can create a Scalar by delivering a value to the constructor. The value type delivered is preserved for the scalar.
An empty constructor prepares a scalar representing zero. Delivering another Scalar copies the value, preserving the value type.

.. code-block:: c++

   Scalar a( 1 );                   // a scalar of IndexType (int)
   Scalar b( 2.5f );                // a scalar of type float
   Scalar c( 3.14 );                // a scalar of type double
   Scalar d( Complex( 2.0, 1.5 ) ); // a scalar of type ComplexDouble
   Scalar e();                      // a scalar representing zero
   Scalar f( c );                   // a new scalar with the value and value type of c

Expressions
-----------

All common binary operators '+', '-', '\*', '/' (also as '+=', '-=', '\*=','/=') and the unary operator '-' are supported.
A mix of different template types is allowed.

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

.. code-block:: c++

   c = max( a, b );
   c = min( a, b );
   c = abs( a );
   c = sqrt( a );


Comparation operators
---------------------

Also the relational operators '==', '!=', '<', '>', '<=', '>=' and the rudimental math functions 'max', 'min', 'abs',
'sqrt' are prepared.

.. code-block:: c++

   bool boolean = ( a == b );
   boolean = ( a != b );
   boolean = ( a < b );
   boolean = ( a > b );
   boolean = ( a <= b );
   boolean = ( a >= b );

Output operator
---------------

For printing a Scalar value you can use the output operator ('<<') or 'getValue<ValueType>()' to print it on your own,
e.g. in printf(). Using the output operator results in 'Scalar(<your_value>)' while you get the plain value with
getValue<ValueType>.

.. code-block:: c++

   // output operator
   std::cout << "my Scalar is: " << a << std::endl;
   
   // getValue
   printf("my Scalar is: %d\n", a.getValue<int>() );