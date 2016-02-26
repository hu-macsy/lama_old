:orphan:

Solution Task 1 (a)
-------------------

.. literalinclude:: ../../examples/lecture/task1a.cpp 
   :language: c++
   :lines: 0-24,37-38
   :emphasize-lines: 21-23

(1) The filename is given as a command-line argument and is the argument of the SparseMatrix-Constructor.
(2) You can get the number of rows by using the method getNumRows().

:download:`Download complete solution Task 1 (a) <../../../solver/examples/lecture/task1a.cpp>`

Solution Task 1 (b)
-------------------

The following Code is an alternative solution for task 1. Explanations for each
line are listened below.

.. literalinclude:: ../../examples/lecture/task1b.cpp 
   :language: c++
   :lines: 0-34,47-48
   :emphasize-lines: 16,18,33

(1) Creation of SparseAssemblyStorage of type double with size 4x4.
(2) Setting some values by using set(). You should only set Non-Zero-Values.
(3) Creation of CSRSparseMatrix of type double and committing the SparseAssemblyStorage.

Setting the right hand side and the solution vector
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../examples/lecture/task1b.cpp 
   :language: c++
   :emphasize-lines: 35,36,37,42,44

(1) Creation of DenseVector rhs of type double and default-values 0.0.
(2) Creation of HostWriteAccess of type double for DenseVector rhs. The Constructor requires a LAMA-Array. You can get it by calling the getLocalValues()-method of your DenseVector.
(3) Setting values of rhs by yourself. The overloaded operator[] makes it easy to handle it.
(4) Release of HostWriteAccesses. Instead of releasing the HostWriteAccess you can use a block { /\* set() here \*/ }. The release()-method will be automatically called of the Destructor at the end of this block.
(5) Creation of DenseVector solution. Default-value is 0.0.

:download:`Download complete solution Task 1 (b) <../../../solver/examples/lecture/task1b.cpp>`

.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_1`", ":doc:`../lecture`", ":doc:`task_2`"
