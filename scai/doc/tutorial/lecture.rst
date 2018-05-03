Lecture
=======

This lecture gives an insight into LAMA and its concepts. It guides you through the LAMA-architecture
and the most important classes, and offers tasks for practical experience. The
practical parts of this lecture illustrate step by step how to create an equation
system, how to solve it and how to run it on different devices. The concepts of 
logging and tracing in LAMA are introduced as well.

Before you start, you should follow the installation
to obtain a working LAMA. With task 0 this lecture gives a step by step introduction 
to the basic usage of LAMA. The following tasks are only giving hints which
functions and classes you should use to complete them. Please consult the
`online API documentation`_ on how to use these classes.

.. _online API documentation: https://libama.org/doxygen/index.html

Requirements
------------

Here we assume that LAMA has already been downloaded, configured, compiled, and
installed on your machine.
The installation directory (**INSTALL_PREFIX** used for the configuration) 
should be set via the environment variable **SCAI_ROOT**.

.. code-block:: bash

   export SCAI_ROOT=<path/to/lama/installation/directory>

The installation directory should contain online documentation:

.. code-block:: bash

   firefox ${SCAI_ROOT}/share/doc/system/html/index.html

The compilation of your LAMA application, e.g. 'simple.cpp' is usually done as
follows:

.. code-block:: bash

   g++ -std=c++11 -I${SCAI_ROOT}/include -o simple simple.cpp -L${SCAI_ROOT}/lib -lscai_lama

For running the executable, it is necessary to include the lib directory of LAMA
into your library path.

.. code-block:: bash

   export LD_LIBRARY_PATH=$SCAI_ROOT/lib:$LD_LIBRARY_PATH

The following example program can be used to verify that compilation, linking 
and running works fine.

:download:`solver <solver.cpp>`

An example ``Makefile`` can be found together with the ``simple.cpp`` example
and the solutions for the lecture tasks at ``${SCAI_ROOT}/doc/user_guide/cpp_source/tutorial/``.

.. code-block:: bash

    ls ${SCAI_ROOT}/share/examples/scai-lama-*/tutorial

.. csv-table:: 
   :header: "Tasks", "Solutions"

   ":doc:`Task 0 <lecture/task_0>`", ""
   ":doc:`Task 1 <lecture/task_1>`", ":doc:`Solution 1 <lecture/solution_task_1>`"
   ":doc:`Task 2 <lecture/task_2>`", ":doc:`Solution 2 <lecture/solution_task_2>`"
   ":doc:`Task 3 <lecture/task_3>`", ":doc:`Solution 3 <lecture/solution_task_3>`"
   ":doc:`Task 4 <lecture/task_4>`", ":doc:`Solution 4 <lecture/solution_task_4>`"
   ":doc:`Task 5 <lecture/task_5>`", ":doc:`Solution 5 <lecture/solution_task_5>`"
   ":doc:`Task 6 <lecture/task_6>`", ":doc:`Solution 6 <lecture/solution_task_6>`"
