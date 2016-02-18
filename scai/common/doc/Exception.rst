Exception
---------

Within LAMA one base class is provided for all LAMA relevant exceptions.

.. code-block:: c++

  class Exception
  {
      ...
  };

This exception class is the base class for all exceptions defined in any SCAI project.
One major feature is that this exception prints the call stack from the place where
the exception has been thrown. This is very helpful to identify the reason for it.

This macros is helpful to throw a general exception:

.. code-block:: c++

   COMMON_THROWEXCEPTION( msg )

This macro can be used to throw a specific exception of a derived class.

.. code-block:: c++

   SCAI_THROWEXCEPTION( exceptionClass, msg )

