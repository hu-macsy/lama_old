Distributions
=============

rowwise --> halo, local

Distribution Typs
-----------------

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
