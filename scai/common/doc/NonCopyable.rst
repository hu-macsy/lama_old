.. _NonCopyable:

NonCopyable
===========

NonCopyable is a class where the default copy operator and the default assignment operator are private.
If a class derives from this class, the default copy operator and the default assignment operator cannot
be used as they rely on the corresponing operators of the base class that are disabled.

.. code-block:: c++

    #include <scai/common/NonCopyable.hpp>

    class MyClass : scai::common::NonCopyable
    {
        ....
    };

    MyClass a, b;
    b = a;            // generates a compile time error
    MyClass c( a );   // generates a compile time error

Default copy and assignment operators should not be used when the class contains member variables
where the ownership might cause problems and where the destructor is not trivial. E.g. a pointer variable will only be copied as a reference and not
deeply so the copied object and the original object contain references to the same object.

It is still possible to provide an own implementation of the copy constructor with its own semantic.

.. code-block:: c++

    class MyClass : scai::common::NonCopyable
    {
        ....

        MyClass( const MyClass& other ) 
        {
            ....
        }
    };

