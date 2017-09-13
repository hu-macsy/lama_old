.. _main-page_benchmark:

#########
Benchmark
#########

***********
Description
***********

The package benchmark provides a framework for running benchmarks in a standardized
way.

Usage
-----

Benchmarks have to be defined a standardized way. They will register themselves
in a registry so that they can be run afterwards in a standardized way with
driver programs and driver scripts.

*********
Reference
*********

.. toctree::

   benchmarkFramework
   benchmarkRunnerFile
   createBenchmark
   inputSet
   publishedBenchmarks

*******
Example
*******

.. code-block:: c++

    lamaSolver stencil3D.frm b.mtx 0 x.mtx

************
Dependencies
************

The framework itself uses the common and logging library.

Python is used to run benchmarks.

************
Related Work
************

None.
