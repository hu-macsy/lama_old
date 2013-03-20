.. _tutorial_solution_task6:

Solution of Task 6
==================

Here is the solution of task 6.

In total there are six logging-levels: TRACE, DEBUG, INFO, WARN, ERROR, FATAL. Where TRACE is the most verbose level.

Setting the environment variable LAMA_LOG with a configuration file:

.. code-block:: bash

   export LAMA_LOG=config

Possible configuration in this configuration file:

.. code-block:: bash

   <root> = WARN
   Solver.IterativeSolver.CG = INFO
   Matrix.CSRSparseMatrix = DEBUG
   Distribution.BlockDistribution = INFO

   
.. csv-table::
   :header: "back to this Task", "Index", "next Task"
   :widths: 330, 340, 330

   ":doc:`task_6`", ":doc:`index`", ":doc:`task_7`"
   