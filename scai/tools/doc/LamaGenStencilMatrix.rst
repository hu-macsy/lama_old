.. _LamaGenStencilMatrix:

********************
lamaGenStencilMatrix
********************

This executable allows to generate a stencil matrix in any file format.

Usage:

.. code-block:: c++

   lamaGenStencilMatrix <filename> <dim> <stencilType> <dimX> [ <dimY> [ <dimZ> ] ]
     filename : name of the output file for matrix, vector
       filename = <id>.mtx -> generates matrix market format, <id>_v.mtx for vector
       filename = <id>     -> generates binary format, <id>.frm for matrix, <id>.frv for vector
       %s in filename is replaced with stencil values, e.g. 2D5P_100_100
     dim = 1, 2, 3  is dimension of stencil
     stencilType = 3 (for dim = 1) 
     stencilType = 5, 9 (for dim = 2) 
     stencilType = 7, 19, 27 (for dim = 3) 

Other options:

.. code-block:: c++

      --SCAI_IO_TYPE_DATA=float|double|ComplexFloat|ComplexDouble to force output in other format
      --SCAI_IO_TYPE_INDEX=int|long to fource output of indexes in other format
      --SCAI_IO_BINARY=0|1 to force formatted or binary format


