.. _Function:

Function/Bind
=============

The C++11 standard provides with the ``function`` a new family of class templates that is similiar
to the function pointer but more general.

Both, function and the traditional function pointer have in common that the signature of the function
is needed for its declaration.

.. code-block:: c++

  int f ( double x, int p )
  {
     ...
  }

  ...

  int ( *foo1 ) ( double, int ) = f;             // traditional function pointer
  std::function<int( double, int )> foo2 = f;    // function wrapper class

The bind routine allows for a function to bind one or more arguments with
a certain argument. If an argument is bound, it does no longer appear in the signature.

.. code-block:: c++

  using namespace std::placeholders;
  std::function<int( double )> foo2a = std::bind( f, _1, 5 );
  std::cout << "Call foo2a( 3.15143 ) = " << foo2a( 3.15143 )  << std::endl;

  std::function<int( int )> foo2b = std::bind( f, 3.15143, _1 );
  std::cout << "Call foo2b( 1 ) = " << foo2b( 1 )  << std::endl;

Placeholders like ``_1`` and ``_2`` are used to define the order of the unbound arguments.
The can also be used to switch arguments.

.. code-block:: c++

  std::function<int( int, double )> foo3 = std::bind( f, _2, _1 );
  std::cout << "Call foo3( 2, 3.15143 ) = " << foo3( 2, 315143 ) << std::endl;

The SCAI libraries use the concept of function and bind at many places. 
Therefore the corresponding C++11 features are used.
