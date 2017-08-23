.. _Partitioning

Partitioning
============

Graph partitioning is a well studied problem in combinatorial scientific computing. 
An important application is the distribution of rows of a sparse matrix for a
matrix vector multiplication where the goals are to balance the load and to minimize 
communication.

LAMA itself does not implement own graph partitioning algorithms but provides
interfaces to well established software packages like Metis.
