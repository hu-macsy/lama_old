OpenMP
------

OpenMP is exploited for the parallelization on one node, either for multiple cores or multiple cpus.

.. code-block:: c++

    #include <omp.h>

    #include <OpenMP.h>

Use of OpenMP.hpp provides dummies for OpenMP routines:
omp_get_num_threads, omp_get_thread_num, omp_get_max_threads
Furthermore: inline routines for atomic operations (ATOMIC is not supported for Complex types).

