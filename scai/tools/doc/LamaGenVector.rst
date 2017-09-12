.. _LamaGenVector:

*************
lamaGenVector
*************

This executable allows to generate vector data any file format.

Usage:

.. code-block:: c++

   lamaGenVector [--SCAI_var=val] outfile_name <size> <val> [matrix_filename]
     outfile_name is filename for vector output
      --SCAI_IO_BINARY=0|1 force formatted or binary output
      --SCAI_TYPE=float|double|LongDouble|ComplexFloat|ComplexDouble|ComplexLong value type
     size is the number of elements in the vector
     val is the value for each entry
      -random each entry is multiplied with a random value from 0..1
     matrix_filename : if set, compute vector as rhs of matrix * vector


