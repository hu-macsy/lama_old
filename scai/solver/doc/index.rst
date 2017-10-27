.. _main-page_solver:

######
Solver
######

***********
Description
***********

The SCAI solver library implements with LAMA's text book syntax linear iterative solvers. The solver palette includes splitting methods like Jacobi or Richardson, krylow subspace methods like CG or GMRES and a simple multigrid method. General stopping criteria are available and preconditioning is preserved.

*****************
Library Reference
*****************

A full list of all prepared solvers can be found :doc:`here <SolverList>`.

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

In the following we explain how to use a solver and solver logging as well as how to set and combine stopping criteria and how to add preconditioning.

=======================    ==========================================
How to                     Description
=======================    ==========================================
:doc:`SolverUsage`         How to use a solver
:doc:`SolverLogger`        How to use to solver logger
:doc:`StoppingCriteria`    How to use stopping criteria
:doc:`Preconditioning`     How to use a preconditioning
=======================    ==========================================

.. :doc:`WritingSolver`    ..   How to write your own solver

.. toctree::
   :hidden:
   
   SolverUsage
   SolverLogger
   StoppingCriteria
   Preconditioning

*******
Example
*******

The following C++ program shows how to use a linear equation solver like the Conjugate Gradient Method.
The input matrix is read from a file, the rhs vector is set to 1 and the initial solution to 0.
:ref:`solver-logging` is defined even it it is optional. 
:ref:`stopping-criteria` is either if a maximal number of iterations is reached or if the L2 norm of the 
residual is smaller than a given value.
The solve method requires the two vectors solution and rhs of the
equation system. Internally the CG calls its iterate method until the defined stopping criterion
is satisfied.

.. code-block:: c++

   #include <scai/lama.hpp>
   #include <scai/solver.hpp>

   using namespace scai::lama;
   using namespace scai::solver;

   int main()
   {
      typedef double ValueType;

      const IndexType maxIter = 20;
      const double eps        = 1e-5;

      CSRSparseMatrix<ValueType> csrMatrix( "gr_30_30.mtx" );

      DenseVector<ValueType> rhs( csrMatrix.getNumRows(), 1.0 );
      DenseVector<ValueType> solution( csrMatrix.getNumColumns(), 0.0 );

      LoggerPtr logger( new CommonLogger ( "myLogger: ",
                                           LogLevel::convergenceHistory,
                                           LoggerWriteBehaviour::toConsoleOnly ) );

      NormPtr norm( new L2Norm() );

      CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
      CriterionPtr it( new IterationCount( maxIter ) );
      CriterionPtr both( new Criterion ( it, rt, Criterion::OR ) );

      CG solver( "simpleExampleCG" );

      solver.setLogger( logger );
      solver.setStoppingCriterion( both );

      solver.initialize( csrMatrix );
      solver.solve( solution, rhs );

      solution.writeToFile( "resultVec.mtx", File::MATRIX_MARKET );
   }

The solver library provides a more complex example program where nearly all parameters
can be set by command line arguments and/or environment variables. 

See :ref:`solver-example`.

.. toctree::
   :hidden:

   SolverExample

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
