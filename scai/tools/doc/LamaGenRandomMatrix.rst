.. _LamaGenRandomMatrix:

*******************
lamaGenRandomMatrix
*******************

This executable allows to generate random matrices in any file format.

Usage:

.. code-block:: c++

    lamaGenRandomMatrix <filename> nrows ncols [ fillrate ]

* filename is is the namoe of the output file for the generated matrix
* The size of the matrix will be nrows x ncols
* The fillrate must be value between 0 and 1 adns specifies the fill rate of the matrix; the value 1 implies
  that the generated matrix is dense.

* `--SCAI_IO_BINARY=<bool>` can be used to force binary or formated output of the matrix data.
  But keep in mind that many matrix file formats only support one kind.

.. code-block:: c++

    lamaGenRandomMatrix matrix_bin.frm 100 100 0.1 --SCAI_IO_BINARY=1
    lamaGenRandomMatrix matrix_form.frm 100 100 0.1 --SCAI_IO_BINARY=0

* `--SCAI_TYPE=<type>` specifies the value type of the matrix.

.. code-block:: c++

    lamaGenRandomMatrix matrix.mtx 5 5 0.2 --SCAI_TYPE=double

    %%MatrixMarket matrix coordinate real general
    5 6 6
    1 6 0.329554488258021
    3 6 -0.213937753079401
    4 1 0.51422645851432
    4 2 -0.608353509202929
    4 3 0.198111211134224
    4 4 0.782382395847125
    
.. code-block:: c++

    lamaGenRandomMatrix matrix.mtx 5 5 0.2 --SCAI_TYPE=ComplexFloat

    %%MatrixMarket matrix coordinate complex general
    5 6 4
    1 6 0.329554 -0.536459
    3 5 -0.213938 0.967399
    4 1 -0.608353 0.686642
    4 3 0.782382 -0.997849


