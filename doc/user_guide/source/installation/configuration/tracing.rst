.. _installation_tracing:

Tracing with LAMA
=================

Within LAMA, many routines have been instrumented at source code level by so-called regions. If tracing is enabled, an
internal subroutine is called at each entry and exit of such a region. These internal subroutines are used for timing
of the regions and/or for the generation of Vampir trace files that can be visualized by Vampir.

By default, tracing is disabled (OFF) and the regions are not instrumented. Tracing can be enabled by the cmake
variable ``LAMA_TRACE_LEVEL`` for timing (TIME) or Vampir tracing (VT).

.. code-block:: bash 

   cmake -DLAMA_TRACE_LEVEL=<OFF | TIME | VT>

Tracing for Time Measurements
-----------------------------

The value ``TIME`` will instrument LAMA in such a way that all regions in the code
are timed.

.. code-block:: bash 

   cmake -DLAMA_TRACE_LEVEL=TIME ...

Beside the source code instrumentation for the installation, take into account that
the timing must be explicitly enabled at runtime by:

.. code-block:: bash 

   export LAMA_TRACE=TIME

Source code instrumentation might add some overhead at runtime due to the calls of
the internal subroutines. When tracing is disabled at runtime, this overhead should
be very low and negligible. Regions in LAMA are usually at a higher lever and never
within time-consuming loops.

Using Vampir/VampirTrace in LAMA
--------------------------------

VampirTrace is an open source library that allows detailed logging of program execution for parallel applications 
using message passing (MPI) and threads (OpenMP, Pthreads). Besides these typical parallelization paradigms, 
VampirTrace is capable of tracing GPU accelerated applications and generates exact time stamps for all GPU related events.

The LAMA library can be instrumented in such a way that at runtime log data is written 
by VampirTrace in the Open Trace Format (OTF), which can be analyzed and visualized by the visualization tool `Vampir`_.
Using VampirTrace is an optional feature of LAMA and must be enabled explicitly.

During the build, source code instrumentation of LAMA for VampirTrace must be enabled in the following way:

.. code-block:: bash 

   export VT_ROOT=<VampirTrace/installation/directory>
   cmake -DLAMA_TRACE_LEVEL=VT ...

The environment variable ``VT_ROOT`` should be set to the directory where the VampirTrace installation is available.

Even if LAMA has been instrumented for VampirTrace, trace files will not be generated automatically at runtime.
Therefore, the overhead of using an instrumented library remains very low. The generation of tracefile must
be enabled explicitly at runtime.

.. code-block:: bash 

    export LAMA_TRACE=vt

Now for a running LAMA application tracefiles will be generated.

.. _Vampir: http://www.vampir.eu

