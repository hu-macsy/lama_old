Matrix Storage
==============

Position of the Diagonal Element
--------------------------------

For some operations on a sparse matrix it is crucial to know the position of the diagonal element in a row. A possible
solution would be to store the index of the diagonal in a separate array. This array is however not necessary if we
follow two rules for all sparse matrix storage schemes

- A diagonal element of a matrix is an element with i = j, where i is the row index and j is the column index. It does
  not matter if the matrix is square or not.
  
- We always store a diagonal element, even if it is zero.

- The diagonal element is stored first in each row, the following elements have an arbitrary order.

These rules do lead to implications if the matrix is distributed. Because a distributed matrix consists of two local
matrices, the local and the halo matrix. A second issue is the fact that the column distribution of a distributed
matrix might not be the same as the row distribution. Because a direct access to the diagonal element is mainly
important for square matrices which have the same distribution for row and columns, we only apply the diagonal first
rule if that is the case. In this cases the diagonal elements are all stored in the local part and therefor the rule
is only applied to the local part of the matrix. If the diagonal first rule is applied to a matrix storage this is
indicated by a boolean flag. Because the halo part of a matrix quite often has a lot of only zero rows, the diagonal
first property would introduce a big overhead in this cases.

.. Data Locality with OpenMP ( First Touch )
.. -----------------------------------------
..
.. Hint Array for zero rows
.. ------------------------

Sparse Matrix Formats
---------------------

This example matrix will show how the different storage formats are handled (the main diagonal will always be
highlighted by **bold** print):

======= ======= ======= =======
**6.0** 0.0     0.0     4.0
7.0     **0.0** 0.0     0.0
0.0     0.0     **9.0** 4.0
2.0     5.0     0.0     **3.0**
2.0     0.0     0.0     1.0
0.0     0.0     0.0     0.0
0.0     1.0     0.0     2.0
======= ======= ======= =======

Compressed Sparse Row Storage Format (CSR)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the CSR format all non-zero-elements are shifted to the left. If the option for shifting the main diagonal to the
front is set to true, the main diagonal elements will stand in the very beginning of each row. Zeroes that are off
the main diagonal are ignored in any case. The CSR saves the number of the saved elements (non-zero-elements plus the
number of zeroes on the main diagonal, if the option is set) (numValues), the number of rows of the matrix (numRows),
the number of columns in the original non-compressed format (numColumns) and arrays for all the non-zero-values and
zeroes in the main diagonal (values), the associated columns for these elements (ja) and the indices of the new
beginnings of a row in *values* as well as the value of *numValues* at the very end (ia).

The CSR format with diagonal element shifting for the example matrix looks like this:

======= ======= =======
**6.0** 4.0     *
**0.0** 7.0     *
**9.0** 4.0     *
**3.0** 2.0     5.0
2.0     1.0     *
*       *       *
1.0     2.0     *
======= ======= =======

.. code-block:: c++

    numValues  = 13
    numRows    =  6
    numColumns =  4
    values     = {  6.0,  4.0,  0.0,  7.0,  9.0,  4.0,  3.0,  2.0,  5.0,  2.0,  1.0,        1.0,  2.0 }
    ja         = {    0,    3,    1,    0,    2,    3,    3,    0,    1,    0,    3,          1,    3 }
    ia         = {    0,          2,          4,          6,                9,         11,   11,         13}

The CSR format without diagonal element shifting looks like this:

======= ======= =======
**6.0** 4.0     *
7.0     *       *
**9.0** 4.0     *
2.0     5.0     **3.0**
2.0     1.0     *
*       *       *
1.0     2.0     *
======= ======= =======

.. code-block:: c++

    numValues  = 12
    numRows    =  7
    numColumns =  4
    values     = {  6.0,  4.0,  7.0,  9.0,  4.0,  2.0,  5.0,  3.0,  2.0,  1.0,        1.0,  2.0 }
    ja         = {    0,    3,    0,    2,    3,    0,    1,    3,    0,    3,          1,    3 }
    ia         = {    0,          2,    3,          5,                8,         10,   10,         12}

The C++ code for the matrix vector multiplication using a CSR matrix:

.. code-block:: c++

    for ( IndexType i = 0; i < numRows; ++i )
    {
        result[i] = 0.0;
        for ( IndexType jj = ia[i]; jj < ia[i+1]; ++jj )
        {
            IndexType j = ja[jj];
            result[i] += values[jj] * v[j];
        }
    }


ELLPACK Storage Format (ELL)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ELL format is much like the CSR format, though the compressed matrix is filled with zeroes to obtain a "shortened"
version of the original. The ELL matrix saves the number of rows it (and the original matrix) have (numRows), the
number of columns the ELL format has which is equal to the length of its longest rows (numValuesPerRow), the original
number of columns (numColumns). The total number of values, including the zeroes that are used as a filler, can be
calculated with *numRows* * *numValuesPerRow*. Additionally the ELL format saves three arrays as well, one for
all the values in the ELL matrix, which are stored in column major order (values), one for the number of non-zero
values (plus the main diagonal zeroes, if needed) (ia) and one for their associated columns (ja). In case of the
"filler"-zeroes, ja points at the last element's column of this row.

The ELL format is used for matrices with about equivalent numbers of non-zero-values in each row.

The ELL format with diagonal element shifting for the example matrix looks like this:

======= ======= =======
**6.0** 4.0     0.0
**0.0** 7.0     0.0
**9.0** 4.0     0.0
**3.0** 2.0     5.0
2.0     1.0     0.0
0.0     0.0     0.0
1.0     2.0     0.0
======= ======= =======

.. code-block:: c++

    numRows         =  7
    numColumns      =  4
    numValuesPerRow =  3
    values     = {  6.0,  0.0,  9.0,  3.0,  2.0,  0.0,  1.0,  4.0,  7.0,  4.0,  2.0,  1.0,  0.0,  2.0,  0.0,  0.0,  0.0,  5.0,  0.0,  0.0,  0.0 }
    ja         = {    0,    1,    2,    3,    0,    0,    1,    3,    0,    3,    0,    3,    0,    3,    3,    1,    3,    1,    3,    0,    3 }
    ia         = {    2,    2,    2,    3,    2,    0,    2 }

The ELL format without diagonal element shifting for the example matrix looks like this:

======= ======= =======
**6.0** 4.0     0.0
7.0     0.0     0.0
**9.0** 4.0     0.0
2.0     5.0     **3.0**
2.0     1.0     0.0
0.0     0.0     0.0
1.0     2.0     0.0
======= ======= =======

.. code-block:: c++

    numRows         =  7
    numColumns      =  4
    numValuesPerRow =  3
    values     = {  6.0,  7.0,  9.0,  2.0,  2.0,  0.0,  1.0,  4.0,  0.0,  4.0,  5.0,  1.0,  0.0,  2.0,  0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0 }
    ja         = {    0,    0,    2,    0,    0,    0,    1,    3,    0,    3,    1,    3,    0,    3,    3,    0,    3,    3,    3,    0,    3 }
    ia         = {    2,    1,    2,    3,    2,    0,    2 }

The C++ code for the matrix vector multiplication using an ELL matrix:

.. code-block:: c++

    // ia is not used here, not necessarily needed to, but useful if non zero values per line vary a lot. 
    for ( IndexType i = 0; i < numRows; ++i )
    {
        ValueType tmp = 0.0;
        for ( IndexType jj = 0; jj < numValuesPerRow; ++jj )
        {
            const IndexType pos = i + jj * numRows;
            const IndexType j = ja[pos];
            tmp += values[pos] * v[j];
        }
        result[i] = tmp;
    }


Jagged Diagonal Storage Format (JDS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JDS format shifts all the non-zero-elements to the very left as the already mentioned formats do. But afterwards
it sorts the rows by length, so the longest row stands on top of the matrix and the shortest at the bottom. Like the
ELL matrix the elements in the *values*-array are entered in column major order. The JDS comes with the integer
*numValues*, *numRows*, *numColumns* and the number of jagged diagonals (which is equal to the number of
columns in the "jagged" Matrix): *numDiagonals*. It contains 5 arrays: One array for the length of each column
(dlg) and one for the length of each row (ilg), the permutation array which shows, where the lines were supposed to
be before the assorting (perm) and the arrays for the elements (values) and their original column indices (ja) as in
the CSR and ELL format.

The JDS format with diagonal element shifting for the example matrix looks like this:

======= ======= =======
**3.0** 2.0      5.0
**6.0** 4.0     *
**0.0** 7.0     *
**9.0** 4.0     *
2.0     1.0     *
1.0     2.0     *
*       *       *
======= ======= =======

.. code-block:: c++

    numValues    = 13
    numRows      =  7
    numColumns   =  4
    numDiagonals =  3
    values     = {  3.0,  6.0,  0.0,  9.0,  2.0,  1.0,  2.0,  4.0,  7.0,  4.0,  1.0,  2.0,  5.0 }
    ja         = {    3,    0,    1,    2,    0,    1,    0,    3,    0,    3,    3,    3,    1 }
    ilg        = { 3, 2, 2, 2, 2, 2, 0 }
    perm       = { 3, 0, 1, 2, 4, 6, 5 }
    dlg        = { 6, 6, 1 }

The JDS format without diagonal element shifting for the example matrix looks like this:

======= ======= =======
2.0      5.0    **3.0**
**6.0** 4.0     *
**9.0** 4.0     *
2.0     1.0     *
1.0     2.0     *
7.0     *       *
*       *       *
======= ======= =======

.. code-block:: c++

    numValues    = 12
    numRows      =  7
    numColumns   =  4
    numDiagonals =  3
    values     = {  2.0,  6.0,  9.0,  2.0,  1.0,  7.0,  5.0,  4.0,  4.0,  1.0,  2.0,  3.0 }
    ja         = {    0,    0,    2,    0,    1,    0,    1,    3,    3,    3,    3,    3 }
    ilg        = { 3, 2, 2, 2, 2, 1, 0 }
    perm       = { 3, 0, 2, 4, 6, 1, 5 }
    dlg        = { 6, 5, 1 }

The C++ code for the matrix vector multiplication using a JDS matrix:

.. code-block:: c++

    for ( IndexType i = 0; i < numRows; i++ )
    {
        ValueType value = 0.0;
        IndexType offset = i;
        for ( IndexType jj = 0; jj < ilg[i]; jj++ )
        {
            value += values[offset] * v[ja[offset]];
            offset += dlg[jj];
        }
        result[perm[i]] = value;
    }

The array ilg is employed when constructing the JDS array. After sorting the rows the array can be easily recomputed
as follows:

.. code-block:: c++

    for ( IndexType i = 0; i < numRows; i++ )
    {
        numValuesInRow = 0;
        for (k = 0; k < numDiagonals; k++)
        {
            if (dlg[k] < i) break;
            ++numValuesInRow;
        }
        ilg[i] = numValuesInRow;
    }


The C++ code for the matrix vector multiplication without the array ilg would look like this:

.. code-block:: c++

    for ( IndexType i = 0; i < numRows; i++ )
    {
        ValueType value = 0.0;
        IndexType offset = i;
        for ( IndexType k = 0; k < numDiagonals; k++ )
        {
            if (dlg[k] < i) break;
            value += values[offset] * v[ja[offset]];
            offset += dlg[k];
        }
        result[perm[i]] = value;
    }

Diagonal Storage Format (DIA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The DIA format extremely differs from the previous ones. It keeps the matrix in order, but extends every diagonal of
the matrix, that contains a non-zero-element. The other diagonals are ignored completely. The extension is done by
adding zeroes "outside" the matrix until all diagonals have the same specific length (numElementsPerDiagonal). This
specific length is either the number of rows (numRows) or the number of columns (numColumns) and depends on which of
these two integer holds the larger value. The number of diagonals is saved as well (numDiagonals) and the total number
of values (including the added zeroes, excluding the zeroes that were "deleted" in order to ignore the unnecessary
diagonals) is calculated like this: *numValues* = *numDiagonals* * *numElementsPerDiagonal*. All the elements
are stored in diagonal major order in an array (values) and another array shows the offset of the main diagonal
(offset). Negative values in the offset array represent diagonals "below" the main diagonal (its original position),
positive values represent diagonals "above" or "right" from the main diagonal.

The DIA format with diagonal element shifting for the example matrix looks like this:

======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= 
**6.0** 0.0     0.0     0.0     0.0     0.0     0.0     *       4.0     *       *       *       *       *       *    
*       **0.0** 0.0     0.0     0.0     0.0     7.0     0.0     *       0.0     *       *       *       *       *   
*       *       **9.0** 0.0     0.0     0.0     0.0     0.0     4.0     *       0.0     *       *       *       *   
*       *       *       **3.0** 0.0     0.0     2.0     5.0     0.0     0.0     *       0.0     *       *       *   
*       *       *       *       0.0     0.0     2.0     0.0     0.0     1.0     0.0     *       0.0     *       *   
*       *       *       *       *       0.0     0.0     0.0     0.0     0.0     0.0     0.0     *       0.0     *   
*       *       *       *       *       *       0.0     1.0     0.0     2.0     0.0     0.0     0.0     *       0.0
======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= 

.. code-block:: c++

    numValues              = 56
    numRows                =  7
    numColumns             =  4
    numDiagonals           =  8
    numElementsPerDiagonal =  7
    values     = {  6.0,  0.0,  9.0,  3.0,  0.0,  0.0,  0.0,
                    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,
                    0.0,  0.0,  0.0,  0.0,  2.0,  0.0,  0.0,
                    0.0,  0.0,  0.0,  2.0,  0.0,  0.0,  2.0,
                    0.0,  0.0,  0.0,  5.0,  0.0,  0.0,  0.0,
                    0.0,  7.0,  0.0,  0.0,  1.0,  0.0,  0.0,
                    0.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,
                    4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 }
    offset     = { 0, -5, -4, -3, -2, -1,  1,  3 }

The DIA format without diagonal element shifting for the example matrix looks like this:

======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= 
0.0     0.0     0.0     0.0     0.0     **6.0** 0.0     *       4.0     *       *       *       *       *       *   
*       0.0     0.0     0.0     0.0     7.0     **0.0** 0.0     *       0.0     *       *       *       *       *   
*       *       0.0     0.0     0.0     0.0     0.0     **9.0** 4.0     *       0.0     *       *       *       *   
*       *       *        0.0     0.0    2.0     5.0     0.0     **3.0** 0.0     *       0.0     *       *       *   
*       *       *       *        0.0    2.0     0.0     0.0     1.0     0.0     0.0     *       0.0     *       *   
*       *       *       *       *       0.0     0.0     0.0     0.0     0.0     0.0     0.0     *       0.0     *   
*       *       *       *       *       *       1.0     0.0     2.0     0.0     0.0     0.0     0.0     *       0.0
======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= ======= 

.. code-block:: c++

    numValues              = 56
    numRows                =  7
    numColumns             =  4
    numDiagonals           =  8
    numElementsPerDiagonal =  7
    values     = {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,
                    0.0,  0.0,  0.0,  0.0,  2.0,  0.0,  0.0,
                    0.0,  0.0,  0.0,  2.0,  0.0,  0.0,  2.0,
                    0.0,  0.0,  0.0,  5.0,  0.0,  0.0,  0.0,
                    0.0,  7.0,  0.0,  0.0,  1.0,  0.0,  0.0,
                    6.0,  0.0,  9.0,  3.0,  0.0,  0.0,  0.0,
                    0.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,
                    4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 }
    offset     = { -5, -4, -3, -2, -1,  0,  1,  3 }

The C++ code for the matrix vector multiplication using a DIA matrix:

.. code-block:: c++

    for ( IndexType i = 0; i < nnu; i++ )
    {
        ValueType accu = 0.0;
        for ( IndexType ii = 0; ii < nd; ++ii )
        {
            const IndexType j = i + offset[ii];
            if ( j >= nnc )
                break;
            if ( j >= 0 )
                accu += data[ i * nd + ii ] * v[j];
        }
        result[i] = accu;
    }

Coordinate Storage Format (COO)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The COO format is a very simple format. In one array, all the row indices for each value are stored (ia), in another
all the column indices (ja) and a third array saves the desired values (values). Logically the arrays all have a
length of *numValues*. 

If the diagonal element shifting is activated, all main diagonal elements are shifted to the beginning of the arrays.
The other elements might be sorted row-wise or column-wise to optimize the access to the values of one row or one
column.

The COO format with diagonal element shifting for the example matrix looks like this:

======= ======= ======= =======
**6.0** **0.0** **9.0** **3.0**
*       *       *       4.0
7.0     *       *       *  
*       *       *       4.0
2.0     5.0     *       *  
2.0     *       *       1.0
*       *       *       *  
*       1.0     *       2.0
======= ======= ======= =======

.. code-block:: c++

    numValues  = 13
    numRows    =  7
    numColumns =  4
    values     = {  6.0,  0.0,  9.0,  3.0,  4.0,  7.0,  4.0,  2.0,  5.0,  2.0,  1.0,  1.0,  2.0 }
    ja         = {    0,    1,    2,    3,    3,    0,    3,    0,    1,    0,    3,    1,    3 }
    ia         = {    0,    1,    2,    3,    0,    1,    2,    3,    3,    4,    4,    6,    6 }

The COO format without diagonal element shifting for the example matrix looks like this:

======= ======= ======= =======
**6.0** *       *       4.0
7.0     *       *       *  
*       *       **9.0** 4.0
2.0     5.0     *       **3.0**
2.0     *       *       1.0
*       *       *       *  
*       1.0     *       2.0
======= ======= ======= =======

.. code-block:: c++

    numValues  = 12
    numRows    =  7
    numColumns =  4
    values     = {  6.0,  4.0,  7.0,  9.0,  4.0,  2.0,  5.0,  3.0,  2.0,  1.0,  1.0,  2.0 }
    ja         = {    0,    3,    0,    2,    3,    0,    1,    3,    0,    3,    1,    3 }
    ia         = {    0,    0,    1,    2,    2,    3,    3,    3,    4,    4,    6,    6 }

The C++ code for the matrix vector multiplication using a COO matrix:

.. code-block:: c++

    for ( IndexType i = 0; i < numRows; ++i )
    {
        result[i] = 0.0;
    }
    for (IndexType k = 0; k < numValues; ++k)
    {
        result[ia[k]] += values[k] * v[ja[k]];
    }

Sparse Matrix Converters
------------------------

In some cases it is necessary to convert the matrix storage formats into other ones. Therefore every storage can be 
converted to CSR and can be initialized from CSR. Through this mechanism every sparse matrix format can be converted
in another one. 
