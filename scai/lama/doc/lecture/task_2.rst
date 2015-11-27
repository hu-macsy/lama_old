:orphan:

Task 2: Implement a CG Iterative Solver
=======================================

The aim of this task is to introduce the linear algebra functionality of LAMA.
Through operator overloading LAMA supports a text book (or pseudocode) like
syntax for basic linear algebra task like multiplication and addition for
scalars, vectors and matrices. Before you can start you need to include the
following headers to have the needed operators available:

::

    #include <scai/lama/expression/MatrixVectorExpressions.hpp>
    #include <scai/lama/expression/VectorExpressions.hpp>

As mentioned before the problem of this task is to write a CG-Solver by creating
a self-provided implementation. A pseudocode implementation of a CG-Solver can
be found `here http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithm`.
As an input you can take all objects created in task 1.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_1`", ":doc:`solution_task_2`", ":doc:`task_3`"