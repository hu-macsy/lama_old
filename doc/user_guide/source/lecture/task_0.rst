:orphan:

Task 0: Solving a system of linear equations
============================================

The aim of this task is to build an application that reads a matrix from file,
initializes a right hand side and solve this equation system with the LAMA
provided CG solver.

First we need include the headers for the types and operations we are going to
use later. **lama.hpp** needs to be included always make sure it is included
first.

::

   // Always include lama.hpp first
   #include <lama.hpp>

   #include <lama/matrix/CSRSparseMatrix.hpp>
   #include <lama/DenseVector.hpp>
   
   // includes operators (+,*) for book syntax
   #include <lama/expression/all.hpp>

   #include <lama/solver/CG.hpp>
   #include <lama/solver/criteria/ResidualThreshold.hpp>
   #include <lama/norm/L2Norm.hpp>

   #include <iostream>

   using namespace lama;

Next we are reading the input matrix from the file provided via a command line
argument. The CSRSparseMatrix provides a constructor to read from a matrix
market file for this.

::

   int main( int argc, char* argv[] )
   {
       if ( argc < 2 )
       {
           std::cerr << "No input file specified" << std::endl;
           exit( -1 );
        }

        //Read a sparse matrix from the passed input file
        CSRSparseMatrix<double> m( argv[1] );

Next we are creating our right hand side and our solution vector of appropriate
size. To verify the correctness of the solution we will compute later we are
calculating a right hand side from the exact solution. LAMA uses operator
overloading and proxy objects to support and easy to use and simple interface
to BLAS functionality, that we are using to calculate the vector **rhs**.

::

        std::cout << "Read matrix m : " << m << std::endl;
        IndexType size = m.getNumRows();

        //Create rhs vector
        DenseVector<double> rhs( size , 0.0 );
        std::cout << "Vector rhs : " << rhs << std::endl;

        //Create solution vector
        DenseVector<double> solution( size , 1.0 );
        std::cout << "Vector solution : " << solution << std::endl;

        //Compute the rhs that fits our solution to be able to calculate the error later
        rhs = m*solution;

        //Reset soution to zero so that there is something to solve
        solution = 0.0;

Now we have created a system of linear equations. To solve this we want to use
the CG solver provided by LAMA. To use this we need to create an object of type
CG. The constructor of CG takes a string to identify that solver. This string is
used distinguish different solvers, especially if the state of the solver is
logged, e.g. the norm of intermediate residues. Because CG is an iterative
solver we need to define some kind of stopping criterion. We want our solver to
stop after the residue of the current solution falls below **1E-8**. For this
we create an object of Type ResidualThreshold. The constructor of
ResidualThreshold takes three arguments:

1. The norm that should be applied to the residue
 
2. The threshold to check for
 
3. If the threshold should check of an absolute or a relative residual reduction.
 
After we have created the ResidualThreshold object we need to pass it to our
solver object cgSolver.

::

        //Create a CG solver
        CG cgSolver( "CGTestSolver" );
        
        //Create a stopping criterion for the iterative solver cgSolver
        NormPtr norm = NormPtr ( new L2Norm ( ) );
    	CriterionPtr criterion ( new ResidualThreshold ( norm, 1E-8, ResidualThreshold::Absolute ) );
        cgSolver.setStoppingCriterion( criterion );

Finally we need to initialize our solver with the system matrix and solve the
system.

::

        //Initialize the solver
        cgSolver.initialize( m );
        
        //Solve m * solution = rhs
        cgSolver.solve( solution, rhs );


Because we know the exact solution we can calculate the error and the L2-Norm of
the error.

::

        //calculate the error and its L2-Norm
        DenseVector<double> error( size, 1.0 );

        error = error - solution;

        std::cout << "L2-Norm of error is " << l2Norm( error ) << std::endl;

        return 0;
        }

The solution of task 0 and an example input matrix can be found in
**$LAMA_ROOT/share/tutorial/solutions**.

The source code is also available for download `here`__ 

__ http://libama.sourceforge.net/tutorial/solutions/task0.cpp

and the example input matrix can be downloaded `here`__.

__ http://libama.sourceforge.net/tutorial/solutions/gr_30_30.mtx

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   "-", "-", ":doc:`task_1`"