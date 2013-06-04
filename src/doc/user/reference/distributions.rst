Distributions
=============

LAMA is provided to work on distributed systems from PC clusters to supercomputers. Communication between the processes
is handled by a given library LAMA is build with ( to now: a given MPI implementation, e.g. openMPI, mvapich, ...; a
PGAS backend is in progress). Data management for the communication is operated internally. 

Data distribution is done line-by-line. So one process always holds a full row of a matrix. Additionally a matrix has a
column distribution which divide the partial matrix of one process in a **local** and **halo** part. Regarding the
matrix-vector-multiplication with a vector having the column distribution of the matrix, the local part of the matrix
can be processed without communication of the vector parts on other processes, while the halo part can not be processed
before communication.
Internally these two parts are stored autonomous in two storages, so the calculation can be executed independently and
the communication can be executed asynchronously to the calculation on the local part. 

Distribution Types
------------------

Block Distribution
^^^^^^^^^^^^^^^^^^

.. image:: fileadmin/LAMA/json/_images/blockweise.png

Cyclic Distribution
^^^^^^^^^^^^^^^^^^^

.. image:: fileadmin/LAMA/json/_images/cyclic.png

Redistribute
------------

Expression Rules:

::

    A = B + C

versus

::  

    A( B + C )

Matrix versus Solver Distribution
---------------------------------

Force distribution solver related 

::

    CG.setDistribution( A.getDistributionPtr() )
