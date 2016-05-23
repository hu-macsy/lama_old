.. _main-page_lama:

####
LAMA
####

***********
Description
***********

LAMA stands for Library of Accelerated Math Applications. It realizes an easy-to-use text-book-syntax for writing algorithms using mathematical notitation. It is intended for the fast development of heterogeneous software in the field of linear algebra. By using the three basic data structures **Scalar**, **Vector** and **Matrix** for the formulation of an algorithm shifts the decision of the executing target architecture, distribution strategy, (sparse) storage format and data type to a later time. Implementing algorithms once and choosen the right configuration when implementing different use cases saves the user time and reduce the maintenance of similar codes. Code is only written once and can be executed on any system, reaching from embedded devices to supercomputers fitting the problem.

This flexibilty is enabled by the underlying layers. Due to the library hmemo the memory of LAMA's data structures is managed hidden from the user on any supported target architecture - for choosing one set a execution **Context**. The library dmemo offers distribution strategies for multi-node installations and prepares communication routines addressed by LAMA, so no explicit formulation for data exchanged is needed - for handling multiple nodes set a **Distribution** with an **Communicator**. Besides the Math Kernel libraries facilitate kernel functions for dense and sparse matrix-vector operations on every back-end - for using a specific data type and storage format use one.

*****************
Library Reference
*****************

Basic data structures
---------------------

For developing your own algorithms with LAMA's mathematical notation, you need to know these three basic data structures:

=====================    ==========================================
Class                    Description
=====================    ==========================================
:ref:`scalar`            A single value
:ref:`vector`            A (distributed) vector of values (on different target architectures)
:ref:`matrix`            A (distributed) matrix of values (on different target architectures)
=====================    ==========================================

.. toctree::
   :hidden:
   
   scalar
   vector
   matrix


Usage
-----

For implementing a dedicated application with a choosen data type and storage format, target architecture,and distribution strategy, you need to know how to use and set them.

======================    ==========================================
Class                     Description
======================    ==========================================
:ref:`setValueType`       How to use a value type
:ref:`setStorage`         How to use/set a matrix storage
:ref:`setContext`         How to set an execution context
:ref:`setDistribution`    How to set a distribution with an communicator
======================    ==========================================

.. toctree::
   :hidden:
   
   setContext
   setDistribution
   setStorage
   setValueType

.. remove

.. toctree::
	:hidden:

	expressions
	io

.. *********
.. Relations
.. *********

********
Examples
********

.. *****************
.. Runtime Variables
.. *****************

************
Dependencies
************

LAMA is dependant of all underlying libraries dmemo, sparsekernel, utilskernel, blaskernel, kregistry, hmemo, tasking, tracing, logging and common as well as there external libraries MPI [#f1]_, GPI [#f1]_, MKL [#f1]_/BLAS, cuBlas [#f1]_, cuSparse [#f1]_, OpenMP [#f1]_, CUDA [#f1]_, Java [#f1]_, pThreads, Boost [#f2]_, dl 

.. rubric:: Footnotes

.. [#f1] optional component
.. [#f2] only needed when the used compiler is not capable of needed C++ 11 features

************
Related Work
************

There are a couple of other frameworks that are working on the same field of
interest. Just to name a few of them:
   
`Honei`_
    Hardware oriented numerics, efficiently implemented

`MTL4`_
    Matrix Template Library 4

`ViennaCL`_
    OpenCL-based

`PetSC`_
    Portable, Extensible Toolkit for Scientific Computation

.. _Honei : http://dribbroc.github.com/HONEI/
.. _MTL4 : http://www.simunova.com/de/node/65
.. _ViennaCL : http://viennacl.sourceforge.net/
.. _PetSC : http://www.mcs.anl.gov/petsc/
