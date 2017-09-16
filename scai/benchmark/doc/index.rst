.. _main-page_benchmark:

#########
Benchmark
#########

***********
Description
***********

The benchmark package provides a framework to run benchmarks in a standardized way. It includes a library that
helps to setup benchmarks and to load and register them dynamically. Furthermore, scripts are provided to
run them and to present the results in a convenient way.

Usage
-----

Benchmarks have to be defined a standardized way. They will register themselves
in a registry so that they can be run afterwards in a standardized way with
driver programs and driver scripts.

The script ``BenchmarkRunner`` located in the bin directory of the installation
directory ``${SCAI_ROOT}`` can be used to list all available benchmarks.

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

    BenchMarkrunner -B

===========================  =========================================================
  Benchmark ID               Benchmark Arguments
===========================  =========================================================
Convert_float                ( TargetStorageFormat, SourceStorageFormat )
PMatrixTimesVector           ( MatrixFormat, ValueType )
===========================  =========================================================

.. code-block:: c++

    BenchMarkrunner -I

===========================  =========================================================
    InputSet ID              InputSet Arguments
===========================  =========================================================
File                         ( filename )         
Poisson                      ( dim, npoints, n1[, n2[, n3]] )
Random                       ( size, fillingGrade )
===========================  =========================================================

************
Dependencies
************

The benchmark framework itself uses the common and logging library. Python as script
language is used to run benchmarks.

Other packages might use the benchmark framework to build benchmarks releated to
their functionality.

************
Related Work
************

None.
