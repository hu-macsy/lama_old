.. _main-page_lama:

####
LAMA
####

***********
Description
***********

LAMA stands for Linear Algebra Math Applications. It realizes an easy-to-use text-book-syntax for writing algorithms using mathematical notation. It is intended for the fast development of heterogeneous software in the field of linear algebra. By using the three basic data structures ``Scalar``, ``Vector`` and ``Matrix`` for the formulation of an algorithm shifts the decision of the executing target architecture, distribution strategy, (sparse) storage format and data type to a later time. Implementing algorithms once and choosen the right configuration when implementing different use cases saves the user time and reduce the maintenance of similar codes. Code is only written once and can be executed on any system, reaching from embedded devices to supercomputers fitting the problem.

This flexibilty is enabled by the underlying layers. Due to the library hmemo the memory of LAMA's data structures is managed hidden from the user on any supported target architecture - for choosing one set a execution ``Context``. The library dmemo offers distribution strategies for multi-node installations and prepares communication routines addressed by LAMA, so no explicit formulation for data exchanged is needed - for handling multiple nodes set a ``Distribution`` with a ``Communicator``. Besides the Math Kernel libraries facilitate kernel functions for dense and sparse matrix-vector operations on every back-end - for using a specific data type and storage format use one.

*****************
Library Reference
*****************

Basic data structures
---------------------

For developing your own algorithms with LAMA's mathematical notation, you need to know these data structures:

======================    ==========================================
Class                     Description
======================    ==========================================
:doc:`DenseVector`        A (distributed) dense vector
:doc:`SparseVector`       A (distributed) sparse vector
:doc:`Vector`             A (distributed) vector of values (on different target architectures)
:doc:`Matrix`             A (distributed) matrix of values (on different target architectures)
======================    ==========================================

More informations about the mathematical notation can be found :doc:`here <Expressions>`.

.. toctree::
   :hidden:
   
   DenseVector
   SparseVector
   Vector
   Matrix
   Expressions
   DistributedIO
   FileIO
   FileFormat

Usage
-----

For implementing a dedicated application with a choosen data type and storage format, target architecture,and distribution strategy, you need to know how to use and set them.

======================    ==========================================
How to                    Description
======================    ==========================================
:doc:`SetValueType`       How to use a value type
:doc:`SetStorage`         How to use/set a ``MatrixStorage``
:doc:`SetContext`         How to set a ``Context``
:doc:`SetDistribution`    How to set a ``Distribution`` with a ``Communicator``
:doc:`FileIO`             Read and write of local vector/matrix data
:doc:`DistributedIO`      Read and write of distributed vector/matrix data
:doc:`FileFormat`         Supported file formats.
======================    ==========================================

.. toctree::
   :hidden:
   
   SetValueType
   SetStorage
   SetContext
   SetDistribution

.. **************
.. Relationsships
.. **************

*******
Example
*******

.. literalinclude:: ../../lama/examples/tutorial/example.cpp 
   :language: c++
   :lines: 29-

.. *****************
.. Runtime Variables
.. *****************

************
Dependencies
************

LAMA is dependant of all underlying libraries:

* :ref:`SCAI Dmemo <scaidmemo:main-page_dmemo>`
* :ref:`SCAI SparseKernel <scaisparsekernel:main-page_sparsekernel>`
* :ref:`SCAI UtilsKernel <scaiutilskernel:main-page_utilskernel>`
* :ref:`SCAI BLASKernel <scaiblaskernel:main-page_blaskernel>`
* :ref:`SCAI Kregistry <scaikregistry:main-page_kregistry>`
* :ref:`SCAI Hmemo <scaihmemo:main-page_hmemo>`
* :ref:`SCAI Tasking <scaitasking:main-page_tasking>`
* :ref:`SCAI Tracing <scaitracing:main-page_tracing>`
* :ref:`SCAI Logging <scailogging:main-page_logging>`
* :ref:`SCAI Common <scaicommon:main-page_common>`

And therefore inherit all there external dependencies:

* |Boost|: for additional C++ features [#f2]_

.. |Boost| raw:: html

   <a href="http://www.boost.org/" target="_blank"> Boost </a>

* |PThread|: for asynchronous tasking

.. |PThread| raw:: html
   
   <a href="http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/pthread.h.html" target="_blank"> Pthread </a>

* |Java|: for tasking GUI

.. |Java| raw:: html
   
   <a href="https://www.java.com/en/" target="_blank"> Java </a>

* |CUDA|: for Nvidia GPU programming [#f1]_

.. |CUDA| raw:: html

   <a href="http://www.nvidia.com/object/cuda_home_new.html" target="_blank"> CUDA </a>

* |MKL|: for Intel optimized blas routines [#f1]_

.. |MKL| raw:: html

   <a href="https://software.intel.com/en-us/intel-mkl" target="_blank"> Intel MKL </a>

* |CUBLAS|: for Nvidia optimized blas routines [#f1]_

.. |CUBLAS| raw:: html

   <a href="https://developer.nvidia.com/cublas" target="_blank"> cuBLAS </a>


* |CUSPARSE|: for Nvidia optimized sparse blas routines [#f1]_

.. |CUSPARSE| raw:: html

   <a href="https://developer.nvidia.com/cusparse" target="_blank"> cuSPARSE </a>

* |MPI| (Message Passing Interface): for Interprocess Communication [#f1]_

.. |MPI| raw:: html

	<a href="https://www.mpi-forum.org/docs/docs.html" target="_blank"> MPI </a>

* |Metis|: for Graph Partitioning [#f1]_

.. |Metis| raw:: html

	<a href="http://glaros.dtc.umn.edu/gkhome/views/metis" target="_blank"> Metis </a>

.. rubric:: Footnotes

.. [#f1] optional component
.. [#f2] only needed when the used compiler is not capable of needed C++ 11 features

************
Related Work
************

There are a couple of other frameworks that are working on the same field of interest. Just to name a few of them:
   
|MTL4| (Matrix Template Library 4)

.. |MTL4| raw:: html
	
	<a href="http://http://www.simunova.com/de/node/65" target="_blank"> MTL4 </a>

|ViennaCL|

.. |ViennaCL| raw:: html

	<a href="http://viennacl.sourceforge.net/" target="_blank"> ViennaCL </a>

|Paralution| 

.. |Paralution| raw:: html

   <a href="http://www.paralution.com/" target="_blank"> Paralution </a>

|PETSc| (Portable, Extensible Toolkit for Scientific Computation)

.. |PETSc| raw:: html

	<a href="http://www.mcs.anl.gov/petsc/" target="_blank"> PETSc </a>
