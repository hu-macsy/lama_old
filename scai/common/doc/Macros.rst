.. _Macros:

``SCAI_COMMON_COUNT_NARGS``
---------------------------

The macro ``SCAI_COMMON_COUNT_NARG`` can be used to get the number of arguments.

.. code-block:: c++

   #include <scai/common/macros/count.hpp>


   #define HOST_TYPES int, float, double, long double
   ...
   std::cout << SCAI_COMMON_COUNT_NARG( HOST_TYPES ) << " types are supported: " << HOST_TYES << std::endl;

This is the output:


.. code-block:: c++

   4 types are supported: int, float, double, long double


``SCAI_UNUSED``
---------------

Most C++ compilers give warnings when an argument of a function is never
used within the function. And this is indeed helpful as this is mostly
a bad coding style.

Nevertheless, there might be situations, when it cannot be avoided. Within
LAMA, the signature of a kernel function is the same for all devices on which
the routine is implemented and there it might be the case that one argument
is not needed for a certain implementation on some device. 

.. code-block:: c++

   #include <scai/common/macros/unused.hpp>

   void func( float& a, const int SCAI_UNUSED( n ), const int m )
   {
       a = static_cast<float>( m );
   }

``SCAI_SYSTEM_CALL``
--------------------

Many system functions return a status value that indicates whether the operation 
has been successful or not.

.. code-block:: c++

   int rc = pthread_init( ... );
   if ( rc != 0 )
   {
      COMMON_THROWEXCEPTION( ... )
   }

A corresponding macro can be used that will also contain the line number of the 
system call. 

.. code-block:: c++

   #include <scai/common/system_call.hpp>
   
   SCAI_SYSTEM_CALL( pthread_init( ... ), msg )
   SCAI_SYSTEM_CALL_NOTHROW( pthread_init( ... ), msg )

The macro ``SCAI_SYSTEM_CALL_NOTHROW`` is similiar but gives only an
error message and does not throw an exception. It might be used in 
implementations of destructor methods where it is bad practice to throw an 
exception.

``SCAI_UNSUPPORTED``
--------------------

Not every operation within LAMA might be available on each platform for all kind
of precisions, storage formats or distributions among the heterogeneous nodes of the platform. 
Generally speaking, two strategies are possible:

* an exception is thrown for an unsupported feature,
* or the operation is done less efficient by using workarounds like 
  redistributions or other devices
  that might involve conversion and communication overhead.

Within LAMA, the second strategy is followed wherever possible. Either the operation itself might
not be cruicial for the execution time, or the operation might be implemented more efficiently 
by corresponding extensions of the methods later. At all these places, the ``SCAI_UNSUPPORTED``
macro is used. It guarantees that these inefficient operations can be used during development and
can be identified for optimization of the application.

.. code-block:: c++

   #include <scai/common/macros/unsupported.hpp>
   
   SCAI_UNSUPPORTED( msg )

The action of this macro depends on the status of the environment variable
``SCAI_UNSUPPORTED``:

* ``IGNORE`` nothing
* ``WARN`` prints a warning with the messagem
* ``ERROR`` throws an exception with msg as error message

``SCAI_COMMON_UNIQUE_NAME``
---------------------------

The following macro can be used to generate an unique name by joining 
a prefix, the line number, and a postfix.

.. code-block:: c++

   #include <scai/common/macros/unique_name.hpp>

   int SCAI_COMMON_UNIQUE_NAME( lock, var );  // lock7var 
   int SCAI_COMMON_UNIQUE_NAME( lock, var );  // lock8var

This macro itself is used within other macros that require a temporary variable.

``SCAI_COMMON_LOOP``
--------------------

The macro ``SCAI_COMMON_LOOP`` is very useful to apply one macro to multiple arguments.
Therefore it avoids a great deal of code repitition.

.. code-block:: c++

   #include <scai/common/macros/loop.hpp>

   #define MY_MACRO( x ) std::cout << x << std::endl;

   MY_MACRO( 1 )
   MY_MACRO( 2 )
   MY_MACRO( 3 )

   // here is the shorter way

   SCAI_COMMON_LOOP( MY_MACRO, 1, 2, 3 )

The macro ``SCAI_COMMON_LOOP_LVL2`` is similiar but has one additional argument that
can be used within the macro over which will be looped.

.. code-block:: c++

   #define MY_MACRO1( OUT, x ) OUT << x << std::endl;
   
   SCAI_COMMON_LOOP_LVL2( std::err, MY_MACRO1, 1, 2, 3, 2, 1 )
   SCAI_COMMON_LOOP_LVL2( std::cout, MY_MACRO1, 5, 6, 7 )

The macro ``SCAI_COMMON_LOOP_LVL3`` is also availabe when two additional
arguments are needed.

.. code-block:: c++

   #define MY_MACRO2( ARG1, ARG2, x ) ARG1 << ARG2 << x << std::endl;
   
   SCAI_COMMON_LOOP_LVL3( std::err, 1, MY_MACRO2, 1, 2, 3, 2, 1 )
   SCAI_COMMON_LOOP_LVL3( std::cout, 2, MY_MACRO2, 5, 6, 7 )

The combination of the macros is very helpful for nesting.

.. literalinclude:: ../examples/DemoLoopMacro.cpp
   :language: c++
   :lines: 30-

