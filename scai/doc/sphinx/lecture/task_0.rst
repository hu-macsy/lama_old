:orphan:

Task 0: Solving a system of linear equations
--------------------------------------------

The aim of this task is to build an application that reads a matrix from a file,
initializes a right hand side and solves this equation system with the LAMA
provided CG solver.

First we need to include the headers for the types and operations we are going to
use later. **lama.hpp** needs to be included. Always make sure it is included
first.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 3-18

Next we are reading the input matrix from the file provided via a command line
argument. The CSRSparseMatrix provides a constructor to read from a matrix
market file for this.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 20-29

Then we create our right hand side and our solution vector of appropriate
size. To verify the correctness of the solution we are
calculating a right hand side from the exact solution. LAMA uses operator
overloading and proxy objects to support and easy to use and simple interface
to BLAS functionality, that we are using to calculate the vector **rhs**.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 30-45

Now we have created a system of linear equations. To solve it we want to use
the CG solver provided by LAMA and so we start with creating an object of type
CG. The constructor of CG takes a string to identify that solver. This string is
used to distinguish different solvers, especially if the state of the solver is
logged, e.g. the norm of intermediate residues. Because CG is an iterative
solver, we need to define some kind of stopping criterion. We want our solver to
stop after the residue of the current solution falls below **1E-8**. For this
we create an object of Type ResidualThreshold. The constructor of
ResidualThreshold takes three arguments:

1. The norm that should be applied to the residue
 
2. The threshold to check for
 
3. If the threshold should check of an absolute or a relative residual reduction.
 
After we have created the ResidualThreshold object we need to pass it to our
solver object cgSolver.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 47-52

Finally we need to initialize our solver with the system matrix and solve the
system.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 54-57

Because we know the exact solution we can calculate the error and the L2-Norm of
the error.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 59-67

Task 0 and an example input matrix can be found in
**<project-root>/examples/lecture**.

:download:`Download source file Task 0 <../../../solver/examples/lecture/task0.cpp>`

:download:`Download input matrix <../../../solver/examples/lecture/gr_30_30.mtx>`

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   "-", "-", ":doc:`task_1`"
