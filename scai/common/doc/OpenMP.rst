OpenMP
======

OpenMP is exploited for the parallelization on one node, either for multiple cores or multiple CPUs.
OpenMP is mainly based on macros that are enabled by a certain compiler flag.
If any of the OpenMP API routines is called, the file ``omp.h`` must be included.

.. code-block:: c++

    #include <omp.h>

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int nthreads  = omp_get_num_threads();
        ...
    }

Unfortunately, this include statement makes the use of OpenMP mandatory and the code does not
compile any more if the OpenMP compiler flag is not set.
In order to keep OpenMP optional for the compilation of the SCAI libraries, the file
``OpenMP.hpp`` of the common project should be included.

.. code-block:: c++

    #include <scai/common/OpenMP.hpp>
    ...

The file OpenMP.hpp provides dummies for the corresponding OpenMP routines like
``omp_get_num_threads``, ``omp_get_thread_num``, ``omp_get_max_threads``.
Furthermore it contains some inline routines for atomic operations as ``ATOMIC`` is not supported for Complex types.

