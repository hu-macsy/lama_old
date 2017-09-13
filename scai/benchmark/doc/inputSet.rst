Input sets
==========

.. Introduction

Random
------

The **Random** input set takes one or two parameters:

 - size: The number of rows of the matrix M x M and the vectors that will be created.
 - fillingDegree: The filling degree of the matrix (]0.0, 1.0]); default = 1.0.

This input set creates a matrix and two vectors with randomly distributed values of the value (1 / size).

File
^^^^

The **File** input set takes one parameter.

 - file: The path to a matrix file.

It reads the matrix out of the file, which either has to be in the Matrix Market format (\*.mtx) or in the AMG format (\*.frm).
Furthermore it creates two vectors, the first vector has as much values as the matrix columns. The other one is the product of the first vector and the matrix.

Poisson
^^^^^^^

The **Poisson** input set takes one parameter.

 - type: The type of the matrix.

Type should have the following format: **aDbP_dimX_dimY_dimZ** valid values are:

+--------+-----+------+-------+
| a =    | 1   | 2    | 3     |
+--------+-----+------+-------+
| b =    | 3   | 5, 9 | 7, 27 |
+--------+-----+------+-------+
| dimX = | >=1 | >=1  | >=1   |
+--------+-----+------+-------+
| dimY = | 1   | >=1  | >=1   |
+--------+-----+------+-------+
| dimZ = | 1   | 1    | >=1   |
+--------+-----+------+-------+

LAMASpVdotDVInputSet
^^^^^^^^^^^^^^^^^^^^

The **LAMASpVdotDV** input set takes two or three parameters separated with ",":

 - dimension: The dimension of the vectors.

 - fillingDegree: The filling degree (]0.0, 1.0]).

 - seed: The seed for the random generator.

This Input Set only creates one sparse vector and an dense vector.

**!**This Input Set only creates two vectors but no Matrix not for benchmarking matrix operations **!**
