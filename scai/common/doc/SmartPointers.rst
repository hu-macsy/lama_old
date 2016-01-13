Smart Pointers
--------------

The philosophy of using smart pointers in LAMA is explained in the developers guide.
As smart pointers free the allocated data in their destructor, the code becomes 
more robust as a forgotten free does not matter and when an exception is thrown.

Smart pointers are now provided by C++11. For backward compatibility it is possible
to use smart pointers of the Boost libarary. For abstraction, we have corresponding
type definitions in the common library.

.. code-block:: c++

   #include <common/shared_ptr.hpp>
   #include <common/unique_ptr.hpp>

   common::shared_ptr<T> X ( new T (...) );
   common::unique_ptr<T> Y ( new T (...) );
   common::scoped_array<T> A ( new T [10] );


