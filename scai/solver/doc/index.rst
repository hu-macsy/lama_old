.. _main-page_solver:

######
Solver
######

***********
Description
***********

The SCAI solver library implements with LAMA's text book syntax linear iterative solvers. The solver palette includes splitting methods like Jacobi or Richardson, krylow subspace methods like CG or GMRES and a simple multigrid method. Preconditioning is preserved and general stopping criteria a available.

*****************
Library Reference
*****************

A full list of all prepared solvers can be found :doc:`here <SolverList>`

.. toctree::
   :hidden:
   
   SolverList

.. Basic data structures
.. ---------------------

.. =======================    ==========================================
.. Class                      Description
.. =======================    ==========================================
.. :doc:`Solver`                
.. :doc:`IterativeSolver`
.. =======================    ==========================================

.. toctree::
   :hidden:

   Solver
   IterativeSolver

Usage
-----

In the following we explain how to use a solver in general, solver logging, set and combine stopping criteria and add preconditioning.

=======================    ==========================================
How to                     Description
=======================    ==========================================
:doc:`SolverUsage`         How to use a solver
:doc:`SolverLogger`        How to use to solver logger
:doc:`StoppingCriteria`    How to use stopping criteria
:doc:`Preconditioning`     How to use a preconditioning
.. :doc:`WritingSolver`       How to write your own solver
=======================    ==========================================

.. toctree::
   :hidden:
   
   SolverUsage
   SolverLogger
   StoppingCriteria
   Preconditioning

*******
Example
*******

.. code-block:: c++

   #include <scai/lama.hpp>
   #include <scai/solver.hpp>

   using namespace scai::lama;
   using namespace scai::solver;

   int main()
   {
      typedef double ValueType;

      IndexType size = 30;

      DenseVector<ValueType> rhs( size, 1.0 );
      DenseVector<ValueType> solution( size, 0.0 );
      CSRSparseMatrix<ValueType> csrMatrix( "gr_30_30.mtx" );

      LoggerPtr logger( new CommonLogger ( "myLogger",
                                           LogLevel::convergenceHistory,
                                           LoggerWriteBehaviour::toConsoleOnly ) );

	   CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
	   CriterionPtr it( new IterationCount( lamaconf.getMaxIter() ) );
	   CriterionPtr both( new Criterion ( it, rt, Criterion::OR ) );

	   solver->setLogger( logger );
	   solver->setStoppingCriterion( rt );

	   solver->initialize( csrMatrix );
	   solver->solve( solution, rhs );

      solution.writeToFile( "resultVec.mtx", File::MATRIX_MARKET );
   }

************
Dependencies
************

SCAI solver is dependant of all underlying libraries:

* :ref:`SCAI LAMA <scailama:main-page_lama>`
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

* Compiler supporting Intel MIC Architecture for using the Xeon Phi Coprocessor

* |MKL|: for Intel optimized blas routines [#f1]_

.. |MKL| raw:: html

   <a href="https://software.intel.com/en-us/intel-mkl" target="_blank"> Intel MKL </a>

* |CUBLAS|: for Nvidia optimized blas routines [#f1]_

.. |CUBLAS| raw:: html

   <a href="https://developer.nvidia.com/cublas" target="_blank"> CUBLAS </a>


* |CUSPARSE|: for Nvidia optimized sparse blas routines [#f1]_

.. |CUSPARSE| raw:: html

   <a href="https://developer.nvidia.com/cusparse" target="_blank"> CUSPARSE </a>

* |MPI| (Message Passing Interface): for Interprocess Communication [#f1]_

.. |MPI| raw:: html

	<a href="https://www.mpi-forum.org/docs/docs.html" target="_blank"> MPI </a>

* |GPI| (Global Adress Programming Interface): for Interprocess Communication [#f1]_

.. |GPI| raw:: html

	<a href="http://www.gpi-site.com/gpi2" target="_blank"> GPI </a>

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
