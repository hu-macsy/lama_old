Task 7: Performance Visualization with Vampir (optional)
========================================================

Vampir is a tool to measure and visualize the performance of instrumented programs.
During the build process of LAMA it is possible to switch on source
code instrumentation that makes possible the generation of trace files at runtime
(see ":doc:`../installation/tracing`").
This is supported by corresponding calls of the VampirTrace library.
So by using LAMA you can also generate Vampir trace files of your
application. Especially for MPI parallel runs, the communication between the
processes is traced to identify communication patterns and to get communication
statistics.

By default, tracing is deactivated at runtime. Tracing can be activated by
setting the environment variable ``LAMA_TRACE``.

.. code-block:: bash

   export LAMA_TRACE=vt

Furthermore, you can also enable tracing for additional running threads 
that will be started for asynchronous computations.

.. code-block:: bash

   export LAMA_TRACE=vt:thread

When using CUDA, CUDA runtime and device activities can be recorded by
setting the following environment variable (of VampirTrace):

.. code-block:: bash

   export VT_GPUTRACE=yes

By using VampirTrace you can take advantage of all the features offered
by this trace library. Very convenient for performance analysis is the
hardware counter support. 

.. code-block:: bash

   papi_avail     ! shows supported hardware performance counters
   export VT_METRICS=PAPI_L2_DCR:PAPI_FP_OPS

After a successful run of your application, the generated ``.otf`` files can be
visualized by Vampir.

The most interesting methods of LAMA have been instrumented to identify them as
regions in the trace file. LAMA uses corresponding macros to define code
sections as regions. Each region has a name that should be unique for the whole
application. These macros can also be used in user applications.

Lets try it: With the makro LAMA_REGION("CUDA region") create a new region at
the beginning of your program before setting the CUDA-Context. Use the makros
LAMA_REGION_START("CG-region") and LAMA_REGION_END("CG-region") around the
solve-method() of your CG-Solver to create a second region. Run the program MPI
parallel and open the created \*.otf file with Vampir. To be able to see all of
your regions it is advisable to increase the zoom of the timeline. You can find
a specific region by its name by using Ctrl+F.

.. csv-table:: 
   :header: "previous", "Solution", "next"
   :widths: 330, 340, 330

   ":doc:`task_6`", ":doc:`solution_task_7`", "-"
