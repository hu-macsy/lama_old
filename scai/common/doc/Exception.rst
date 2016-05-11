.. _Exception:

Exception
=========

For all projects one base class is provided that is used as base class for all derived exceptions.

.. code-block:: c++

  class Exception : std::exception
  {
      ...
  };

One major feature is that this exception can print the call stack from the place where
the exception has been thrown. This is very helpful to identify the location in the call
hierarchy where the exception has been thrown.

This exception class is the base class for all exceptions defined in any SCAI project.
Derived exception classes will not have any new functionality but are helpful to act on different
situations. The common subproject itself provides alread the following derived exceptions: 
AssertException, NotSupportedValueTypeException, UnsupportedException.  

.. code-block:: c++

  try 
  {
      ...
  }
  catch ( scai::common::AssertException& ex )
  {
     std::cout << "Assertion failed. " << ex.what() << std::endl;
     return -1;
  }
  catch ( scai::common::UnsupportedValueException& ex )
  {
     std::cout << "Assertion failed. " << ex.what() << std::endl;
     return -1;
  }

The common subproject provides two macros that should be used whenever an exception is 
thrown. The advantage of the macro is that the source code location (filename, line number) will
be included in the message of the exception and that the operator << can be used for the 
composition of messages.

.. code-block:: c++

   #include <scai/common/macros/throw.hpp>
   ...
   COMMON_THROWEXCEPTION( "This is the reason why we will throw an exception, value = " << value )

This macro can be used to throw a specific exception of a derived class.

.. code-block:: c++

   #include <scai/common/macros/throw.hpp>
   ...
   SCAI_THROWEXCEPTION( exceptionClass, msg )

