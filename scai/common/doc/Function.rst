Function/Bind
-------------

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

  int ( *foo1 ) ( double, int ) = f;       // traditional function pointer
  function<int( double, int )> foo2 = f;    // function wrapper class

The bind routine allows for a function to bind one or more arguments with
a certain argument. If an argument is bound, it does no longer appear in the signature.

.. code-block:: c++

  function<int( double )> foo2a = bind( f, _1, 5 );
  std::cout << "Call foo2a( 3.15143 ) = " << foo2a( 3.15143 )  << std::endl;

  function<int( int )> foo2b = bind( f, 3.15143, _1 );
  std::cout << "Call foo2b( 1 ) = " << foo2b( 1 )  << std::endl;

Placeholders like ``_1`` and ``_2`` are used to define the order of the unbound arguments.
The can also be used to switch arguments.

.. code-block:: c++

  function<int( int, double )> foo3 = bind( f, _2, _1 );
  std::cout << "Call foo3( 2, 3.15143 ) = " << foo3( 2, 315143 ) << std::endl;

The SCAI libraries use the concepts of function and bind at many places. As 
the standard C++11 is not yet available on all machines, the correspoding concepts
of the Boost libraries might be used. By using ``function`` and ``bind`` classes
of the namespace ``scai::common``, it is possible to keep the code independent of the
actual solution that is taken on a certain platform.

.. code-block:: c++

  #include <scai/common/function.hpp>
  #include <scai/common/bind.hpp>

  using namespace scai::common;

  function<int( double, int )> foo2 = &f;    // function wrapper class
  function<int( double )> foo2a = bind( f, _1, 5 );


