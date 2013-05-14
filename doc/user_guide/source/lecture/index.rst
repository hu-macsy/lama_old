*******
Lecture
*******

This lecture gives an insight into the Library for Accelerated Math
Applications(LAMA) and its concepts. It guides you through the LAMA-architecture
and the most important classes, and offers tasks for practical experience. The
practical parts of this lecture illustrate step by step how to create an equation
system, how to solve it and how to run it on different devices. The concepts of 
logging and tracing in LAMA are introduced as well.

Before you start, you should follow the :doc:`../installation/index`
to obtain a working LAMA. With task 0 this lecture gives a step by step introduction 
to the basic usage of LAMA. The following tasks are only giving hints which
functions and classes you should use to complete them. Please consult the
`online API documentation`_ on how to use these classes.

.. _online API documentation: http://libama.sourceforge.net/doc/index.html

Requirements
============

Here we assume that LAMA has already been downloaded, configured, compiled, and
installed on your machine.
The installation directory (**INSTALL_PREFIX** used for the configuration) 
should be set via the environment variable **LAMA_ROOT**.

.. code-block:: bash

   export LAMA_ROOT=<path/to/lama/installation/directory>

The installation directory should contain online documentation
that is also available online `here`__.

__ http://libama.sourceforge.net/doc/index.html

.. code-block:: bash

   firefox <project-root>/doc/doxygen/html/index.html

The compilation of your LAMA application, e.g. 'simple.cpp' is usually done as
follows:

.. code-block:: bash

   g++ -I $LAMA_ROOT/include -L $LAMA_ROOT/lib -llama -llog4lama -o simple simple.cpp

For running the executable, it is necessary to include the lib directory of LAMA
into your library path.

.. code-block:: bash

   export LD_LIBRARY_PATH=$LAMA_ROOT/lib:$LD_LIBRARY_PATH

The following example program can be used to verify that compilation, linking
and running works fine.

:doc:`Solver <../tutorial/example_solver>`

An example **Makefile** can be found together with the **simple.cpp** example
and the solutions for the lecture tasks at **<project-root>/doc/user_guide/cpp_source/tutorial/**.

.. H4H Tutorial Remarks
.. ====================

.. To run the tutorial on nova you need to log in to the head node of nova and submit an
.. interactive job to the gpus queue. Please allocate 2 cpus so that all tutorial
.. participants can get free resources and we are able to run MPI parallel jobs
.. later in this tutorial. The tutorial will also need the two modules mentioned
.. below.

.. code-block:bash

   qsub -Iq gpus -lnodes=1:ppn=2
   module load bullxmpi/bullxmpi-1.0.2
   module load intel_compiler/12.0.2.137
   export LAMA_ROOT=/home_nfs/h4h/LAMA/lama

.. csv-table:: 
   :header: "Tasks", "Solutions"

   ":doc:`Task 0 <task_0>`", ""
   ":doc:`Task 1 <task_1>`", ":doc:`Solution 1 <solution_task_1>`"
   ":doc:`Task 2 <task_2>`", ":doc:`Solution 2 <solution_task_2>`"
   ":doc:`Task 3 <task_3>`", ":doc:`Solution 3 <solution_task_3>`"
   ":doc:`Task 4 <task_4>`", ":doc:`Solution 4 <solution_task_4>`"
   ":doc:`Task 5 <task_5>`", ":doc:`Solution 5 <solution_task_5>`"
   ":doc:`Task 6 <task_6>`", ":doc:`Solution 6 <solution_task_6>`"
   ":doc:`Task 7 <task_7>`", ":doc:`Solution 7 <solution_task_7>`"
