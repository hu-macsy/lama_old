:orphan:

Task 0: Solving a system of linear equations
--------------------------------------------

The aim of this task is to build an application that reads a matrix from a file,
initializes a right hand side and solves this equation system with the LAMA
provided CG solver.

First we need to include the headers for the types and operations we are going to
use later. ``lama.hpp`` needs to be included. Always make sure it is included
first. At the beginning of the main program we verify that there is at least one
command line argument.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 36-63

Next we read the input matrix from the file provided via a command line
argument. All Matrix classes provide a constructor that reads 
the matrix data from a file. If only the file name is specified, the suffix
is used to determine which file format is used. Most likely we use files that have the
suffix ``.mtx`` that stands for the Matrix Market file format.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 64-67

For the solver example we create our right hand side and our solution vector of appropriate
size. To verify the correctness of the solution we are
calculating a right hand side from the exact solution. LAMA uses operator
overloading to suppport matrix-vector multiplicatioan that is
used here to construct the vector ``rhs``.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 68-77

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
 
After we have created the ResidualThreshold object we pass it to our
solver object cgSolver.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 78-83

Finally we initialize our solver with the system matrix and solve the
system.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 84-88

Because we know the exact solution we can calculate the error and the L2-Norm of
the error.

.. literalinclude:: ../../../solver/examples/lecture/task0.cpp 
   :language: c++
   :lines: 89-93

Task 0 and an example input matrix can be found in
``${SCAI_ROOT}/share/examples/scai-solver-*/lecture`` or downloaded here.

:download:`Download source file Task 0 <../../../solver/examples/lecture/task0.cpp>`

:download:`Download input matrix <../../../solver/examples/lecture/gr_30_30.mtx>`

**Exercise**: Download a symmetric positive definite matrix
from Matrix Market [http://math.nist.gov/MatrixMarket/searchtool.html] and
use it as input file for the solver program presented here.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   "-", "-", ":doc:`task_1`"
