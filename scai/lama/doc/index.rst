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

For developing your own algorithms with LAMA's mathematical notation, you need to know these three basic data structures:

======================    ==========================================
Class                     Description
======================    ==========================================
:doc:`Scalar`             A single value
:doc:`Vector`             A (distributed) vector of values (on different target architectures)
:doc:`Matrix`             A (distributed) matrix of values (on different target architectures)
======================    ==========================================

More informations about the mathematical notation can be found :doc:`here <Expressions>`.

.. toctree::
   :hidden:
   
   Scalar
   Vector
   Matrix
   Expressions

Usage
-----

For implementing a dedicated application with a choosen data type and storage format, target architecture,and distribution strategy, you need to know how to use and set them.

======================    ==========================================
Class                     Description
======================    ==========================================
:doc:`SetValueType`       How to use a value type
:doc:`SetStorage`         How to use/set a ``MatrixStorage``
:doc:`SetContext`         How to set a ``Context``
:doc:`SetDistribution`    How to set a ``Distribution`` with a ``Communicator``
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

* |MPI| (Message Passing Interface)

.. |MPI| raw:: html

	<a href="https://www.mpi-forum.org/docs/docs.html" target="_blank"> MPI </a>

* |GPI| (Global Adress Programming Interface)

.. |GPI| raw:: html

	<a href="http://www.gpi-site.com/gpi2" target="_blank"> GPI </a>

* |Metis| Graph Partitioning Software

.. |Metis| raw:: html

	<a href="http://glaros.dtc.umn.edu/gkhome/views/metis" target="_blank"> Metis </a>

dmemo, sparsekernel, utilskernel, blaskernel, kregistry, hmemo, tasking, tracing, logging and common as well as there external libraries MPI [#f1]_, GPI [#f1]_, MKL [#f1]_/BLAS, cuBlas [#f1]_, cuSparse [#f1]_, OpenMP [#f1]_, CUDA [#f1]_, Java [#f1]_, pThreads, Boost [#f2]_, dl 

.. rubric:: Footnotes

.. [#f1] optional component
.. [#f2] only needed when the used compiler is not capable of needed C++ 11 features

************
Related Work
************

There are a couple of other frameworks that are working on the same field of interest. Just to name a few of them:
   
|MTL4| Matrix Template Library 4

.. |MTL4| raw:: html
	
	<a href="http://http://www.simunova.com/de/node/65" target="_blank"> MTL4 </a>

|ViennaCL| OpenCL-based

.. |ViennaCL| raw:: html

	<a href="http://viennacl.sourceforge.net/" target="_blank"> ViennaCL </a>

|PETSc| Portable, Extensible Toolkit for Scientific Computation

.. |PETSc| raw:: html

	<a href="http://www.mcs.anl.gov/petsc/" target="_blank"> PETSc </a>
