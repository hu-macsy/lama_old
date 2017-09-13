.. _LamaMatrixConvert:

*****************
lamaMatrixConvert
*****************

This executable allows to convert matrix file from one format to another format.

Usage:

.. code-block:: c++

   lamaMatrixConvert infile_name outfile_name

The file format of the input or output file is given by its suffix, e.g. frm for the SAMG format,
mtx for the Matrix Market format and so on.

*   --SCAI_TYPE=<data_type> is data type of input file and used for internal representation
*   --SCAI_IO_BINARY=0|1 to force formatted or binary output file
*   --SCAI_IO_TYPE_DATA=<data_type> is data type used for file output


