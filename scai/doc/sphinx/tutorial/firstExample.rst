First LAMA Example Program
--------------------------

The following C++ program shows a very simple example program of how to use LAMA.

.. literalinclude:: ../../../lama/examples/tutorial/simple.cpp 
   :language: c++
   :lines: 34-

:download:`Download source file <../../../lama/examples/tutorial/simple.cpp>`

The include file lama.hpp contains some definitions how far assertions, logging and tracing statements
are compiled into your code. The definitions will be the same as used for the installation.

Usually you have to include the class definition file for each LAMA class that you are
using. As we use objects of class DenseVector and Scalar, we have to include the corresponding files.

For more informations about includes of LAMA see :doc:`here <include>`.
