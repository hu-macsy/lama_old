.. _hcdk:

Heterogeneous Computing Development Kit
---------------------------------------

LAMA's base module facilitates three issues:

* A consistent data usage on heterogeneous devices via dedicated read and write accesses within the memory management.
* Decisions about the execution context in the application are separated from implementing the kernels by the kernel management. 
* Asynchronous execution of these kernels, memory transfer and communication is handled by the tasking layer. This combination leads to a clean software development, accomplishing a good maintainability on the user's side with minimal runtime overhead.

* :ref:`SCAI Common - Basic Concepts <scaicommon:main-page_common>`
* :ref:`SCAI Logging - Logging Macros <scailogging:main-page_logging>`
* :ref:`SCAI Tracing - Tracing Macros <scaitracing:main-page_tracing>`
* :ref:`SCAI Tasking - Asynchronous Tasks <scaitasking:main-page_tasking>`
* :ref:`SCAI Hmemo - Heterogeneous Memory Architecture <scaihmemo:main-page_hmemo>`
* :ref:`SCAI Kregistry - Generic Kernel Registry <scaikregistry:main-page_kregistry>`
