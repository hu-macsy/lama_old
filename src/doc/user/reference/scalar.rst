Scalar
======

The class Scalar represents a multi precision scalar. It can be one of these template types: INDEX_TYPE, FLOAT or DOUBLE.

You can create a Scalar by delivering a value to the constructor. An empty constructor prepare a scalar representing zero.

.. code-block:: c++

   Scalar a( 1 );    // a scalar of index type
   Scalar b( 2.5f ); // a scalar of type float 
   Scalar c();       // a scalar representing zero

All common binary operators '+', '-', '*', '/' (also as '+=', '-=', '*=','/=') and the unary operator '-' are supported.
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

Also the relational operators '==', '!=', '<', '>', '<=', '>=' and the rudimental math functions 'max', 'min', 'abs',
'sqrt' are prepared.

.. code-block:: c++

   //relational operators
   bool boolean = ( a == b );
   boolean = ( a != b );
   boolean = ( a < b );
   boolean = ( a > b );
   boolean = ( a <= b );
   boolean = ( a >= b );
   
   // math functions
   c = max( a, b );
   c = min( a, b );
   c = abs( a );
   c = sqrt( a );
   
For printing a Scalar value you can use the output operator ('<<') or 'getValue<ValueType>()' to print it on your own,
e.g. in printf(). Using the output operator results in 'Scalar(<your_value>)' while you get the plain value with
getValue<ValueType>.

.. code-block:: c++

   // output operator
   std::cout << "my Scalar is: " << a << std::endl;
   
   // getValue
   printf("my Scalar is: %d\n", a.getValue<int>() );