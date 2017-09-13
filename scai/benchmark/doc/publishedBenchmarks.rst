Publication Benchmarks
======================

To run the benchmarks presented in the publications about LAMA check out LAMA from our git repository, configure and
build it. How to configure and build LAMA is described in the :doc:`../installation/build`.
Take care that you set **CMAKE_BUILD_TYPE=Release** to obtain an optimized build of LAMA.

For licensing issues we can only provide a binary version of an unoptimized single process AMG setup library `here`__.

__ http://libama.sourceforge.net/benchmark/libamglama.so

It is however possible to use this library for a distributed memory execution of the solution phase of AMG. The fact
that the provided setup library is not optimized does not affect the performance of the solution phase on an matrix
hierarchy created with this library. To allow a distributed memory execution of the solution phase the system matrix is
replicated on every process and each process runs the whole setup. This leads to a quite high memory consumption so it
might not be possible to run all benchmarks that we have presented in our publications depending on the available main
memory.

A file which contains most of the presented benchmarks can be downloaded `here`__.

__ http://libama.sourceforge.net/benchmark/LAMAAMGSolution.beru

To run all benchmarks from this file execute the following:

.. code-block:: bash

   export LAMA_AMG_SETUP_LIBRARY=<path-to-libamglama.so>
   export BENCHMARK_LIBRARY_PATH=<path-to-directory-with-liblama_benchmarks.so>
   BenchmarkRunner -f LAMAAMGSolution.beru

The python script BenchmarkRunner has several options, to list all available options call

.. code-block:: bash

   BenchmarkRunner --help

A description of the file **LAMAAMGSolution.beru** is provided with inline comments. A further introduction to LAMA can
be found in our :doc:`../lecture/index`. If there are any problems in executing the benchmarks please `contact`_ us.

.. _`contact`: mailto:lama@scai.fraunhofer.de
