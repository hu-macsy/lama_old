.. _ReductionOp:

ReductionOp
===========

The following enumeration type specifies the different kind of binary operators
that can be used in reduction operators.

=========  =================================
Name       Operation
=========  =================================
COPY       x = y
ADD        x += y
MULT       x \*= y
MIN        x = min( x, y )
MAX        x = max( x, y )
ABS_MAX    x = max( x, abs(y) )
=========  =================================

The enum class is used in order to have one common function with an addtional op argument instead
of individual functions for each kind of operator.

