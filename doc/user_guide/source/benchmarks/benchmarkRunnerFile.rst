Structure of the Benchmark Runner Files
=======================================

The Benchmark Runner Files are formatted according to the YAML Specification.
The Benchmark Runner File contains 2 Lists:

- BENCHMARKS

- INPUT_SETS

Each List has one or more Elements. The elements are annotated as in the
following example:

.. code-block:: bash

   BENCHMARKS:
    - Benchmark1
    - Benchmark2
   INPUT_SETS:
    - inputSet1(args)
    - inputSet2(args)

As you can see in the example there are two Benchmarks specified and two input
sets. This will result in the following Tests.

- Benchmark1 with inputSet1
- Benchmark1 with inputSet2
- Benchmark2 with inputSet1
- Benchmark2 with inputSet2

Furthermore there are three variables in the file:

- REP_MIN
- TIME_MIN
- REP_CTRL

These variables are annotated as follows:

.. code-block:: bash

   REP_MIN: x
   TIME_MIN: y
   REP_CTRL: z

This Configuratin would result in the benchmarks beeing lauched at least **z**
times to avoid aberrations. One Benchmark will execute its execute() function at
least **x** times. But the benchmark will run at least **y** seconds.

.. code-block:: bash

   TODO LAUNCH_CONF

Input Sets
----------

You can specify arguments, passed to an Input Set, with brackets:

.. code-block:: bash

   INPUT_SETS:
    - inputSet1(400 500)
    - inputSet2(arument)

This will pass the Arguments **400** and **500** to **inputSet1** and
**argument** to **inputSet2** for a list on the specific Input Sets and their
arguments you should have a look at :doc:`inputSet`
