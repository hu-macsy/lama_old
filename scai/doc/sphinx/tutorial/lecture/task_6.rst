:orphan:

Task 6: Enable Logging to track the moving of data to and from the GPU
----------------------------------------------------------------------

LAMA includes sophisticated logging facilities. There are different gradations
of logging: TRACE, DEBUG, INFO, WARN, ERROR, and FATAL. You can regulate the logging
output by setting the environment variable SCAI_LOG with a configuration file.
In a configuration file it is possible to specify logging levels for each
created logger. In general a logger exists for each class with the same
name as the class. The log level of a super class is inherited by derived
classes, to mimic this in the logger configuration the class hierarchy is
encoded in the logger name. E.g. the name of the logger of the class DenseVector
is Vector. DenseVector, because DenseVector implements the Vector interface.
The default log level for all loggers is WARN.

Now its your turn: Set the environment variable with a configure file and change
the logging level of your CSRSparseMatrixLogger (e.g. "Matrix.CSRSparseMatrix =
DEBUG"). If you are using the LAMA Provided CG solver you can activate logging
for it with "Solver.IterativeSolver.CG = INFO".
For a detailed explanation on loggin refer to :ref:`scailogging:main-page_logging`.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_5`", ":doc:`solution_task_6`", "-"