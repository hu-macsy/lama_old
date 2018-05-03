:orphan:

Task 1: Setup a system of linear equations
------------------------------------------

Linear algebra requires working with objects like scalars, vectors and matrices.
LAMA includes classes, which represent those objects and offers the possibility
to handle with different sparse matrix formats ( e.g. CSR, ELL, JDS). This task
introduces some possibilities to setup such matrices and vectors.

Setting up a CSRSparseMatrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

LAMA currently supports the sparse matrix formats CSR, ELLPACK and JDS through
the classes CSRSparseMatrix, ELLSparseMatrix and JDSSparseMatrix. During this
tutorial we will focus on the widely used CSR format.

Lets get started with a good practice: the first essential step is to know how
to create those objects. 

The first possibility and probably the most easiest way to create a matrix is
to read it from file. This possibility has already been shown in the previous
task. 

.. code-block:: c++

    CSRSparseMatrix<ValueType> matrix ( <filename> );

If you want to create a matrix with your own data you should use a 
``MatrixAssemblyAccess`` for this. Such an access offers an efficient
possibility to assemble your data without knowing the quantity of
elements or the order of them. With the release of the access,
the assembled data is converted to the required sparse matrix format.
As we will see later, this approach works also for distributed matrices.
During the access elements are inserted by using the
push method. 

.. code-block:: c++

    CSRSparseMatrix<ValueType> matrix( numRows, numColumns );
    
    MatrixAssemblyAccess mAccess( matrix );
    for ( ... )
    {
        mAccess.push( rowPos, colPos, value );
    }
    mAccess.release();   // now the assembled data is converted to the CSR format

**Exercise**: Write some code that fills the matrix with data.
Take care that the matrix becomes symmetric positive definite so
that it can be used for an CG Iterative Solver works. We propose the
discretization of a Possion Equation.

Setting the right hand side and the solution vector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the first part of this task we set up the matrix. Now we will do the same for (dense) vectors.
For our equation system we need a righthandside-vector and a solution-vector.
Now, your exercise is to create two objects of type DenseVector and fill them
with elements. 

In contrary to a matrix, a dense vector might be initialized
with a value (as shown in the previous task).
To fill a DenseVector with other data, you can use a VectorAssemblyAccess. 

.. code-block:: c++

    DenseVector<ValueType> vector( size, ValueType( 0 ) );
    
    VectorAssemblyAccess vAccess( vector );
    for ( ... )
    {
        vAccess.push( pos, value );
    }
    vAccess.release();   // now the assembled data is filled in the dense vector.

Example Data
^^^^^^^^^^^^

For the exercise you might use the following matrix and rhs data:

.. math::

    matrix = 
  \left(\begin{matrix} 
    1  & 0 & 0 & 0 & 0 & 0 & 0  \\
    -s & (1+2s) & -s & 0 & 0 & 0 & 0 \\
    0  & -s & (1+2s) & -s & 0 & 0 & 0 \\
    0  & 0 & -s & (1+2s) & -s & 0 & 0 \\
    0  & 0 & 0 & -s & (1+2s) & -s & 0 \\
    0  & 0 & 0 & 0 & -s & (1+2s) & -s \\
    0  & 0 & 0 & 0  & 0 & 0 & 1 \\
    \end{matrix}\right) 
    \;
    rhs = 
  \left(\begin{matrix} Tleft \\
    0 \\
    0 \\
    0 \\
    0 \\
    0 \\
    Tright \end{matrix}\right)    

A good choice for the value s is 0.5.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_0`", :doc:`Solution Task 1 <solution_task_1>`, ":doc:`task_2`"
