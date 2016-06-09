.. _lap:

Linear Algebra Package
----------------------

The Linear Algebra Package facilitates the development of (sparse) numerical algorithms for various application domains. Code can be written in text-book-syntax as

.. code-block:: c++

	y = A * x

(where x and y are vectors and A is a ­matrix). Due to the underlying layers, the problem formulation is handled independently of the implementation details regardless of the target architecture and distribution strategy as memory management and communication is processed internally. Furthermore, with load balancing between different components and asynchronous execution, full system performance can be obtained.
In addition, LAMA offers various iterative solvers like Jacobi or CG methods, that can be used directly or preconditioned, with a combination of several user-definable ­stopping criteria. Furthermore, the integration of a custom-built solver is straightforward.

* :ref:`SCAI LAMA - LAMA core <scailama:main-page_lama>`
* :ref:`SCAI Solver - Iterative Solver <scaisolver:main-page_solver>`