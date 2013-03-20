Input Sets
==========

Random
------

The Random Input Set takes one or two parameters:

- size: The number of rows of the matrix M x M and the Vectors that will be created.
- fillingGrade: The fillinggrade of the matrix. The fillinggrade defaults to 1.0 (]0.0, 1.0])

This input set creates a Matrix and two Vectors, with randomly distributed
values of the value (1 / size).

File
^^^^

The File Input Set takes one parameter.

- file: The Path to a Matrix File.

It reads the Matrix out of the File, which either has to be in the Matrix Market
Format (\*.mtx) or in the amg Format (\*.frm).
Furthermore it creates two Vectors, the First Vector has as much Values as the
Matrix columns. The other one is the Product of the first Vector and the Matrix.

Poisson
^^^^^^^

The Poisson Input Set takes one parameter.

- type: The type of the Matrix.

type should have the following format: **aDbP_dimX_dimY_dimZ** valid Values are:

+--------+----+-----+------+
| a =    | 1  | 2   | 3    |
+--------+----+-----+------+
| b =    | 3  | 5, 9| 7, 27|
+--------+----+-----+------+
| dimX = | >=1| >=1 |  >=1 |
+--------+----+-----+------+
| dimY = | 1  | >=1 |  >=1 |
+--------+----+-----+------+
| dimZ = | 1  | 1   |  >=1 |
+--------+----+-----+------+

LAMASpVdotDVInputSet
^^^^^^^^^^^^^^^^^^^^

The LAMASpVdotDV Input Set takes 2 or three Parameters separated with ",":

- dimension: The dimension of the Vectors

- fillingGrade: The fillinggrade (]0.0, 1.0])

- seed: The seed for the random generator

This Input Set only creates one Sparse Vector and an Dense Vector.

**!**This Input Set only creates two vectors but no Matrix not for benchmarking matrix operations **!**