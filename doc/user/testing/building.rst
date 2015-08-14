Building and Executing Unit Tests
=================================

Running the unit tests is very helpful to identify problems at an early stage.
Therefore it is recommended to run the LAMA unit tests after successful compilation.

The unit tests will be built automatically by calling ``make`` in the build directory.
The following executables will be available in the build directory:

* ``lama_test`` in directory ``test`` (will be generated always)
* ``lama_cuda_test`` in directory ``test/cuda`` (will be generated only if CUDA is available and enabled)
* ``lama_dist_test`` in directory ``test/distributed`` (will be generated only if MPI is available and enabled)

In a first step, you should run the unit tests only as a single process (but with
multiple OpenMP threads) and without using any GPU. Here are the commands::

  cd <build-directory>/test
  export SCAI_LOG=ERROR
  export LAMA_UNSUPPORTED=IGNORE
  export LAMA_TEST_CONTEXT=Host
  ./lama_test

* ``SCAI_LOG`` should be set to ``ERROR`` to avoid logging of exceptions that might be thrown 
  by all unit tests that test for failures.
* ``LAMA_UNSUPPORTED`` should be set to ``IGNORE`` to avoid warnings for matrix conversions 
  (here we test functionality and not efficiency).
* ``LAMA_TEST_CONTEXT=Host`` restricts the tests to run only on the CPU.

By default, the tests will use all available processors for the OpenMP parallelized code.
You can set explicitly the number of used processors as follows::

  export OMP_NUM_THREADS=1
  ./lama_test
  export OMP_NUM_THREADS=4
  ./lama_test

If CUDA has been enabled (is default if CUDA is available) and your system has at least one GPU
you can run also the tests on the GPU (you can choose the device by the environment variable
``LAMA_DEVICE``, by default device 0 is taken).

.. code-block:: bash

  export LAMA_TEST_CONTEXT=CUDA
  export LAMA_DEVICE=0
  ./lama_test
  ./cuda/lama_cuda_test

Note: the executable ``lama_cuda_test`` in the subdirectory ``cuda`` contains some special tests
only for the GPU and will not be generated if CUDA is disabled.

Parallel execution of the tests can be started via ``mpirun`` (if MPI has been enabled and found 
during the configuration).

.. code-block:: bash

  export LAMA_TEST_CONTEXT=Host
  mpirun -np <number_of_processors> ./lama_test
  mpirun -np <number_of_processors> ./distributed/lama_dist_test

Note: the executable ``lama_dist_test`` in the subdirectory ``distributed`` contains some
special tests only for MPI and will only be generated if MPI is available and enabled.

If you run the parallel tests on multiple nodes make sure that the environment variables
have the same values on all available nodes. For OpenMPI, this is done as follows:

.. code-block:: bash

  mpirun -x SCAI_LOG=ERROR -x LAMA_UNSUPPORTED=IGNORE
         -x LAMA_TEST_CONTEXT=Host -machine_file=<file> ./lama_test

Parallel execution of the tests on multiple GPUs works also fine as long as one node has
exactly one GPU.

.. code-block:: bash

  mpirun -x SCAI_LOG=ERROR -x LAMA_UNSUPPORTED=IGNORE
         -x LAMA_TEST_CONTEXT=CUDA -x LAMA_DEVICE=0 -machine_file=<file> ./lama_test

If you need a solution where MPI processes require different values for the LAMA_DEVICE
you need a hand-coded script for starting the tests (we will work on a more convenient 
solution).
