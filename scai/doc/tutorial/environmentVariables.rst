:orphan:

Environment variables
---------------------

This is a summary of environment variables that are currently used for LAMA.

SCAI_LOG
^^^^^^^^

The variable ``SCAI_LOG`` specifies the detail of logging or the logger configuration file.

.. code-block:: bash

    export SCAI_LOG=WARN  ( is default )
    export SCAI_LOG=INFO
    export SCAI_LOG=DEBUG
    export SCAI_LOG=<filename>

For specific logging of certain classes a configuration file is required. The entries in this file could be
as follows:

.. code-block:: bash

    Matrix = DEBUG
    Matrix.DenseMatrix = INFO

SCAI_TRACE
^^^^^^^^^^

The variable ``SCAI_TRACE`` can be used to enable the collection of runtime data.

.. code-block:: bash

    export SCAI_TRACE=OFF                   ( is default )
    export SCAI_TRACE=time                  ( collects timing of all instrumented regions )
    export SCAI_TRACE=time:PREFIX=myTiming  ( give trace files a specific prefix )

More detailed information about tracing can be found in the :ref:`SCAI Tracing<scaitracing:main-page_tracing>` module.

SCAI_UNSUPPORTED
^^^^^^^^^^^^^^^^

The variable ``SCAI_UNSUPPORTED`` specifies how the LAMA library deals with operations
that are not efficiently implemented. These might be operations 
that will not be executed at the intended location (e.g. on the CPU instead on the GPU) 
or not in the desired matrix format (e.g. implicit conversion to CSR format and back).

.. code-block:: bash

    export SCAI_UNSUPPORTED=WARN (default)
    export SCAI_UNSUPPORTED=ERROR
    export SCAI_UNSUPPORTED=IGNORE

* ``WARN`` will print always a warning for such an operation
* ``ERROR`` will throw an exception (and so it might terminate if exception is not caught)
* ``IGNORE`` will not give any information or action for such operations.

SCAI_DEVICE
^^^^^^^^^^^

This variable specifies the default device for creation of a CUDA Context in case operations should be
executed on a GPU.

Example

.. code-block:: bash

    export SCAI_DEVICE=1

SCAI_NP
^^^^^^^

This variable specifies the configuration of a processor array. It is used for Poisson input setgenerators
that generates sparse matrices for 2D and 3D problems. 

Examples

.. code-block:: bash

    export SCAI_NP="2 4 2"
    export SCAI_NP=2x2
    export SCAI_NP=2_2_4

The product of the number of processors in each dimension must be the same as the number of processors
on which the application will run.

If the variable is not set, LAMA will take its own factorization of the available processors that fits
best to the given problem size.

Note: Please keep in mind that distributions of matrices are always one-dimensional row distributions and
therefore the environment variable has no influence.
