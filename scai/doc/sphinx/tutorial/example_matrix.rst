Matrix Example
--------------

Multiplication of a CSRSparseMatrix and a DenseVector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we simply multiply a sparse matrix in CSR format and a dense vector.
Our multiplication looks like this:   

.. math::
  \left(\begin{matrix} r_0 \\
    r_1 \\
    r_2 \\
    r_3 \\
    r_4 \\
    r_5 \\
    r_6 \end{matrix}\right) =
  \left(\begin{matrix} 6 & 0  & 0 & 4 \\
    7 & 0 & 0 & 0 \\
    0 & 0 & -9.3 & 4 \\
    2 & 5 & 0 & 3 \\
    2 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 1 & 0 & 2 \end{matrix}\right) *
  \left(\begin{matrix} 6 \\
    4 \\
    7 \\
    -9.3 \end{matrix}\right)

.. literalinclude:: ../../../lama/examples/tutorial/matrix.cpp 
   :language: c++
   :lines: 53-122
    
The result is:

.. math::
  \left(\begin{matrix} r_0 \\
    r_1 \\
    r_2 \\
    r_3 \\
    r_4 \\
    r_5 \\
    r_6 \end{matrix}\right) = 
  \left(\begin{matrix} -1.2 \\
    42 \\
    -102.3 \\
    4.1 \\
    2.7 \\
    0 \\
    -14.6 \end{matrix}\right)    
    
The full example program can be found here :download:`matrix.cpp <../../../lama/examples/tutorial/matrix.cpp>`
	