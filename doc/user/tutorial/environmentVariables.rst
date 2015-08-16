:orphan:

Environment variables
=====================

This is a summary of environment variables that are currently used for LAMA.

SCAI_LOG
--------

The variable ``SCAI_LOG`` specifies the detail of logging or the logger configuration file.

::

    export SCAI_LOG=WARN  ( is default )
    export SCAI_LOG=INFO
    export SCAI_LOG=DEBUG
    export SCAI_LOG=<filename>

For specific logging of certain classes a configuration file is required. The entries in this file could be
as follows:

::

    Matrix = DEBUG
    Matrix.DenseMatrix = INFO

LAMA_UNSUPPORTED
----------------

The variable ``LAMA_UNSUPPORTED`` specifies how the LAMA library deals with operations
that are not efficiently implemented. These might be operations 
that will not be executed at the intended location (e.g. on the CPU instead on the GPU) 
or not in the desired matrix format (e.g. implicit conversion to CSR format and back).

::

    export LAMA_UNSUPPORTED=WARN (default)
    export LAMA_UNSUPPORTED=ERROR
    export LAMA_UNSUPPORTED=IGNORE

* ``WARN`` will print always a warning for such an operation
* ``ERROR`` will throw an exception (and so it might terminate if exception is not caught)
* ``IGNORE`` will not give any information or action for such operations.

SCAI_DEVICE
-----------

This variable specifies the default device for creation of a CUDA Context in case operations should be
executed on a GPU.

Example::

    export SCAI_DEVICE=1

LAMA_TEST_DEVICE
----------------

This variable is only used for the LAMA unit test to restrict execution of the test on a 
specific device. If not set, tests will run on all devices.

::

    export LAMA_TEST_DEVICE=Host
    export LAMA_TEST_DEVICE=CUDA

LAMA_NP
-------

This variable specifies the configuration of a processor array. It is used for Poisson input setgenerators
that generates sparse matrices for 2D and 3D problems. 

Examples::

    export LAMA_NP="2 4 2"
    export LAMA_NP=2x2
    export LAMA_NP=2_2_4

The product of the number of processors in each dimension must be the same as the number of processors
on which the application will run.

If the variable is not set, LAMA will take its own factorization of the available processors that fits
best to the given problem size.

Note: Please keep in mind that distributions of matrices are always one-dimensional row distributions and
therefore the environment variable has no influence.
