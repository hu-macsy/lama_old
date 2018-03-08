.. _main-page_ipbcls:

##################################
Interior-Point Least Square Solver
##################################

***********
Description
***********

This module provides a solver that solves a least square problem with box constraints.

*****************
Library Reference
*****************


*******
Example
*******

The following C++ program shows how to call the ipbcls solver.

.. code-block:: c++

   #include <scai/lama.hpp>
   #include <scai/ipbcls.hpp>

   using namespace scai;
   using namespace lama;
   using namespace ipbcls;

   int main()
   {
       Matrix<doule>& A = 
       DenseVector<double>& b = ...;    // b.getDistribution() == A.getRowDistribution() 
       DenseVector<double>& x = ...;    // x.getDistribution() == A.getColDistribution()
       DenseVector<double>& lb = ....;  // lb.getDistribution() == A.getColDistribution()
       DenseVector<double>& ub = ....;  // ub.getDistribution() == A.getColDistribution()

       ConstrainedLeastSquares<double> lsq( A );

       lsq.setInnerSolverType( InnerSolverType::NewtonStepCG );

       lsq.setObjectiveTolerance( tol );
       lsq.setMaxIter( 500 );

       try
       {
           lsq.solve( x, b, lb, ub );
       }
       catch ( common::Exception& ex )
       {
           std::cout << "Caught exception: " << ex.what() << std::endl;
           std::cout << "Stop execution." << std::endl;
           return 1;
       }
   }

The full example program can be found as ``lsbc.cpp`` in the example directory
of the module.

************
Dependencies
************

SCAI ipblcs is uses for its interface the LAMA data structures for matrices and solver.
Internally it uses some solvers of the solver library.

* :ref:`SCAI LAMA <scailama:main-page_solver`
* :ref:`SCAI LAMA <scailama:main-page_lama>`

************
Related Work
************

To continue.
