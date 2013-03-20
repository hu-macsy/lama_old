Environment variables
=====================

This is a summary of environment variables that are currently used for LAMA.

LOG4LAMA
---------

This variable specifies the detail of logging.

::

	export LOG4LAMA=WARN  ( is default )
	export LOG4LAMA=INFO
	export LOG4LAMA=DEBUG
	export LOG4LAMA=<filename>

For specific logging of certain classes a configuration file is required. The entries in this file could be
as follows:

::

	Matrix = DEBUG
	Matrix.DenseMatrix = INFO

LAMA_DEVICE
------------

This variable specifies the default device for creation of a CUDA Context in case operations should be
executed on a GPU.

Example:

::

	export LAMA_DEVICE=1

NP4LAMA
-------

This variable specifies the configuration of a processor array. It is used for Poisson input setgenerators
that generates sparse matrices for 2D and 3D problems. 

Examples:

::

	export NP4LAMA="2 4 2"
	export NP4LAMA=2x2
	export NP4LAMA=2_2_4

The product of the number of processors in each dimension must be the same as the number of processors
on which the application will run.

If the variable is not set, LAMA will take its own factorization of the available processors that fits
best to the given problem size.

Note: Please keep in mind that distributions of matrices are always one-dimensional row distributions and
therefore the environment variable has no influence.


LAMA_TEST_NP
------------

Number of processes for tests executed with "make TESTNAME_mpi"
