:orphan:

Solution Task 2
---------------

Here is the solution of task 2. The code demonstrate a self-provided CG-Solver. 

.. literalinclude:: ../../../../solver/examples/lecture/task2.cpp 
   :language: c++
   :lines: 37-109
   :emphasize-lines: 52-53,54-55,57-70

(1) Creation of Vectors, that are needed by your implementation.
(2) Creation of Scalars and epsilon.
(3) Here is the self-provided algorithm of a CG-Solver.

:download:`Download complete solution Task 2 <../../../../solver/examples/lecture/task2.cpp>`

An alternative solution for task 2 is the CG-class of LAMA, which contains the same calculation.

.. literalinclude:: ../../../../solver/examples/lecture/task2a.cpp 
   :language: c++
   :lines: 82-95
   
:download:`Download complete alternative solution Task 2 <../../../../solver/examples/lecture/task2a.cpp>`

**Excursion:**

To see results, it is possible to print out the values of vectors. By the way, an explicit function to print out results does not exist.
But it is possible to use a ReadAccess to get reading access to the values.

.. literalinclude:: ../../../../solver/examples/lecture/task2a.cpp 
   :language: c++
   :lines: 87-92

.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_2`", ":doc:`../lecture`", ":doc:`task_3`"
