Shared and Unique Pointers
--------------------------

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

A scoped array has exactly the same purpose as the unique pointer but is used for arrays instead
of single objects. Its implementation uses the new[] and delete[] operator. 

Shared Pointers in LAMA
^^^^^^^^^^^^^^^^^^^^^^^

LAMA uses heayily shared pointers when ownership is not limited to one place. E.g. the
same distribution might be used for multiple vectors and none of the vectors is preferred to take
the ownership.

The ``std::shared_ptr`` is a new C++11 feature and has exactly the semantic as needed for this purpose.
If it is used, the C++11 support must be explicitly enabled by the C++ compiler.

If the C++ compiler does not support C++11 features, the Boost library offers the same functionality
with ``boost::shared_ptr``.

The SCAI common library offers a class ``common::shared_ptr`` that is used in LAMA. This abstraction
allows using either of the both possibilities without any loss of efficiency.

Unique Pointers in LAMA
^^^^^^^^^^^^^^^^^^^^^^^

LAMA uses unique pointers in all situations where ownership is limited to one scope or where
ownership is exactly defined at one place, e.g. as a member variable of an object. It is not
used when ownership changes and not in standard library containers (see 
:ref:`container`),

The ``std::uniqe_ptr`` is a new C++11 feature and has exactly the semantic as needed. 
If it is used, the C++11 support must be explicitly enabled by the C++ compiler.

If the C++ compiler does not support C++11 features, the ``std::auto_ptr`` might be used that
is very similiar but has some limitations and becomes a deprecated feature in newer versions. 

The SCAI common library offers a class ``common::unique_ptr`` that is used in LAMA. This abstraction
allows using either of the both possibilities without any loss of efficiency. Its functionality
is in the same way limited as the ``std::auto_ptr``.

Scoped Arrays in LAMA
^^^^^^^^^^^^^^^^^^^^^

LAMA uses scoped arrays in all situations where lifetime of the array is limited to one scope or where
ownership is exactly defined, e.g. for a member variable of an object.

The ``std::uniqe_ptr`` is a new C++11 feature and has exactly the semantic as needed for LAMA. 
If it is used, the C++11 support must be explicitly enabled by the C++ compiler.

If the C++ compiler does not support C++11 features, the ``boost::scoped_array`` is used.

The SCAI common library offers a class ``common::scoped_array`` that is used in LAMA. This abstraction
allows using either of the both other possibilities without any loss of efficiency.

::

    {
        int np = getSize();
        scai::common::scoped_array<int> counts( new int[np] );
        scai::common::scoped_array<int> displs( new int[np] );

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
to the container. Unfortunately it requires C++11 support and the former solution of std::auto_ptr cannot be
used in containers. So it is also not supported for common::unique_ptr.

::

    std::vector<Matrix> matrices;                          // copy constructor of Matrix used
    std::vector<Matrix*> matrices1;                        // objects are not freed by destructor
    std::vector<common::unique_ptr<Matrix> > matrices2;    // only possible with std::unique_ptr, but not with std::auto_ptr
    std::vector<common::shared_ptr<Matrix> > matrices3;    // recommended, always 

Therefore, the recommeded solution is using shared pointers in containers as it is done in LAMA for all
allocated objects in containers.
