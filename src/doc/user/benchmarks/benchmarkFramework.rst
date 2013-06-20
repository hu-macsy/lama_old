The LAMA Benchmark Framework
============================

Before you start launching Benchmarks, you should have build LAMA as it is described in :doc:`../installation`.

This page describes how to launch existing benchmarks.
For creating entirely new Benchmarks please refer to :doc:`createBenchmark`

Starting your first Benchmark
-----------------------------

Setting up the Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^

At first we need to set the environment variable **BENCHMARK_LIBRARY_PATH**. The BENCHMARK_LIBRARY_PATH is located in

.. code-block:: bash

   <your-lama-build-path>/bench/

so you have to set your environment variable by

.. code-block:: bash

   export BENCHMARK_LIBRARY_PATH=<your-lama-build-path>/bench/

The BenchmarkRunner
^^^^^^^^^^^^^^^^^^^

List all Benchmarks
"""""""""""""""""""

After we set the environment variable we can list all Benchmarks by executing

.. code-block:: bash

   $ framework/BenchmarkRunner -B

Create a Benchmark Runner File
""""""""""""""""""""""""""""""

Normaly it is not necessary to run all Benchmarks. So you can specify a Benchmark Runner File, which specifies, the
benchmarks to run. To create such a file with all benchmarks enabled execute

.. code-block:: bash

   framework/BenchmarkRunner -b bench.beru

Now you can edit this file and delete the unneeded Benchmarks and Inputsets.
For the Structure of the Benchmark Runner File please refer to :doc:`benchmarkRunnerFile`

Run the Benchmarks
""""""""""""""""""

To run the benchmarks specified in a Benchmark Runner File just execute

.. code-block:: bash

   framework/BenchmarkRunner -f bench.beru

other useful options
""""""""""""""""""""

+------------+---------------------------------------------+
| -o <file>  | Write the Benchmark Results to a csv-file   |
+------------+---------------------------------------------+
| --help     | prints usage information to stdout          |
+------------+---------------------------------------------+
| --longhelp | prints detailed usage information to stdout |
+------------+---------------------------------------------+
