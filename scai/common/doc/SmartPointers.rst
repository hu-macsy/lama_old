.. _SmartPointers:

Smart Pointers
==============

The philosophy of using smart pointers in LAMA is explained in the developers guide.
As smart pointers free the allocated data in their destructor, the code becomes 
more robust as a forgotten free does not matter. Especially in case of exceptions
the compiler takes care that the data is freed.

Smart pointers are now provided by C++11 and are used.

.. code-block:: c++

   #include <memory>

   std::shared_ptr<T> X ( new T (...) );
   std::unique_ptr<T> Y ( new T (...) );
   std::unique_ptr<T[]> A ( new T [10] );


If an object has been allocated dynamically, it must be deleted when it is no more used.
In this sense, somebody takes the ownership and has to take care for the deletion and to
make sure that it is no more used somewhere else. This is not always simple and results
often in programming errors. Another problem is that due to exceptions the code for deletion 
might not be executed and so the memory is not freed.

A shared pointer wraps an allocated object in such a way that ownership has not to be 
observed. I.e. when the last shared pointer for the object in memory is destructed, the
object will be deleted.

A unique pointer should be used when ownership is clearly restricted. Very often this is the
case when the lifetime of the object is in one scope. The unique pointer only takes care that 
the object is freed at the end of the scope or in case of exceptions. It is not recommended
when ownership changes, i.e. in pointer assignments or as arguments or results of functions.

Unique pointers might also be used for an array instead of a single object. In this case the
array type must be specified. Its implementation uses the ``delete[]`` operator instead
of ``delete``.

Shared Pointers
^^^^^^^^^^^^^^^

Shared pointers are the best choice when ownership is not limited to one place. E.g. the
same distribution might be used for multiple vectors and none of the vectors is preferred to take
the ownership.

The ``std::shared_ptr`` is a new C++11 feature and has exactly the semantic as needed for this purpose.

Unique Pointers in LAMA
^^^^^^^^^^^^^^^^^^^^^^^

LAMA uses unique pointers in all situations where ownership is limited to one scope or where
ownership is exactly defined at one place, e.g. as a member variable of an object. It is not
used when ownership changes and not in standard library containers (see 
:ref:`container`),

The ``std::uniqe_ptr`` is a new C++11 feature and has exactly the semantic as needed. 
As it is used, the C++11 support must be explicitly enabled by the C++ compiler.

Scoped Arrays in LAMA
^^^^^^^^^^^^^^^^^^^^^

LAMA uses scoped arrays in all situations where lifetime of the array is limited to one scope or where
ownership is exactly defined, e.g. for a member variable of an object.

The ``std::uniqe_ptr`` is a C++11 feature and has exactly the semantic as needed for LAMA. 
Therefore the C++11 support must be explicitly enabled by the C++ compiler.

.. code-block:: c++

    {
        int np = getSize();
        std::unique_ptr<int[]> counts( new int[np] );
        std::unique_ptr<int[]> displs( new int[np] );

        int displacement = 0;

        for( int i = 0; i < np; i++ )
        {
            counts[i] = static_cast<int>( sizes[i] );
            displs[i] = displacement;
            displacement += counts[i];
        }

        MPI_Scatterv( _, counts.get(), displs.get(), ... )

        // delete[] of counts, displs implicitly by destructor
    }


.. _container:

Containers of Allocated Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Allocated objects are often used in containers (e.g. std::vector) instead of usual objects. This is mandatory
if the copy constructor is not supported for the container type. It might be also more efficient if the 
copy constructor is rather expensive. The problem of ownership and deletion of the objects remains the same
as for all pointer objects. 

Using std::unique_ptr is one of the recommended solutions for containers if ownership is clearly restricted
to the container. The use of std::shared_ptr is also possible and has the advantage that ownership might be
taken over by other objects. As it is more flexible, this approach is more often used in LAMA.

::

    std::vector<Matrix> matrices;                       // copy constructor of Matrix used
    std::vector<Matrix*> matrices1;                     // objects are not freed by destructor
    std::vector<std::unique_ptr<Matrix> > matrices2;    // might be possible
    std::vector<std::shared_ptr<Matrix> > matrices3;    // recommended, always 

Therefore, the recommeded solution is using shared pointers in containers as it is done in LAMA for all
allocated objects in containers.
