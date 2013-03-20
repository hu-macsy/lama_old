Task 6: Enable Logging to track the moving of data to an from the GPU
=====================================================================

LAMA includes sophisticated logging facilities. There are different gradations
of logging: DEBUG, INFO, WARN, ERROR and FATAL. You can regulate the logging
output by setting the environment variable LAMA_LOG with a configuration file.
In a configuration file it is possible to specify logging levels for each
created logger. In general there exists a logger for each class with the same
name than the class. The log level of a super class is inherited by derived
classes, to mimic this in the logger configuration the class hierarchy is
encoded in the logger name. E.g. the name of the logger of the class DenseVector
is Vector. DenseVector, because DenseVector implements the Vector interface.
The default log level for all loggers is WARN.

Now its your turn: Set the environment variable with a configure file and change
the logging level of your CSRSparseMatrixLogger (e.g. "Matrix.CSRSparseMatrix =
DEBUG"). If you are using the LAMA Provided CG solver you can activate logging
for it with "Solver.IterativeSolver.CG = INFO".

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_5`", ":doc:`solution_task_6`", ":doc:`task_7`"