FAQs
====

To submit bugs or for getting help please contact us at lama[at]scai.fraunhofer.de.
Please consider your question is not answered below in our FAQs.

General
-------

**What kind of computer do I need to use LAMA?**

LAMA targets systems that range from desktops and small compute clusters to super computing systems.

Thus, special attention is paid to the parallel scalability and the serial performance of LAMA, facilitating an efficient execution on all of these systems. Currently LAMA supports the applications programming interface OpenMP to exploit Multicore Processors, CUDA and OpenCL (in progress) to exploit accelerators like GPUs and MPI to handle distributed memory systems.

Your system does not need to support all these features. Please refer to our :ref:`software prerequisites<requirements>` to find out the required and optional packages for the LAMA installation.

**What kind of license is LAMA released under?**

For a straightforward and hassle-free integration into commercial software, LAMA is available under the MIT License, a free software license originating at the Massachusetts Institute of Technology.

**How can I contribute to the LAMA development?**

You can contact us via lama[at]scai.fraunhofer.de and confirm our contributor agreement. Then you can send us patches or get access to our git.

**Where can I find documentation?**

We provide two kind of documentation: a user documentation and a developer documentation. This user documentation will be build with the library installation with the property BUILD_DOC=ON (default). The System Documentation can be found |SysDoc| or you can build it on your own by calling make doxygendoc in your build directory.

.. |SysDoc| raw:: html

	<a href="https://test.libama.org/doxygen/index.html" target="_blank"> here </a>

Usage
-----

**What kind of matrix input/output formats are supported by LAMA?**

LAMA supports two main kinds of matrix input/output formats yet: matrix market and SAMG format (options: formatted, binary, xdr). For creating a matrix read in from file just give the matrix constructor a string with the path to the file (eg.: matrix.mtx).

**Can I add my own sparse matrix format or solver to LAMA?**

Yes, you can ;-) LAMA can be easily extended with new data structures for sparse matrices if this is needed to achieve optimal performance on special hardware or problems with unique features. Solvers can be written easyly on your own - mostly you have to implement a new ``iterate`` step in LAMA's text-book-syntax.
.. An example for creating a new solver is given :ref:`here<scaisolver:solver_writingSolver>`.

Execution
---------

**How can I debug my LAMA application?**

Compiled in debug mode you can debug your LAMA application with any debugger your familiar with. LAMA also provides an extra logging concept to check sizes, etc. of common variables. How to use logging is described :ref:`here<scailogging:main-page_logging>`.

**How can I analyse the performance of my LAMA application?**

Within LAMA, many routines have been instrumented at source code level by so-called regions. If tracing is enabled, an internal subroutine is called at each entry and exit of such a region. These internal subroutines are used for timing of the regions and/or for the generation of files in OTF (Open Trace Format) that can be visualized e.g. by Vampir. By default, tracing is disabled (OFF) and the regions are not instrumented. Tracing can be enabled by the cmake variable LAMA_TRACE_LEVEL for timing (TIME) or Vampir tracing (VT).
For more details about tracing in LAMA please refer to :ref:`Tracing with LAMA<scaitracing:main-page_tracing>`.

**Does the LAMA logging and tracing affect the efficiency of my application?**

Logging can slow down an application because of massive prints if it is used on very low level. Therefore LAMA provides several logging level and can be switch on on class level so you can decide what and how much information you want to receive. You should use it for debugging and testing your application only either. For measurements you should switch Logging off.

Even if LAMA has been instrumented for VampirTrace, trace files will not be generated automatically at runtime. Therefore, the overhead of using an instrumented library remains very low.