.. _BinaryOp:

BinaryOp
===========

The following enumeration type specifies the different kind of binary operators
that can be used in set and reduction operators.

=========  =================================
Name       Operation
=========  =================================
COPY       x = y
ADD        x += y
SUB        x -= y
MULT       x \*= y
DIVIDE     x /= y
MIN        x = min( x, y )
MAX        x = max( x, y )
ABS_MAX    x = max( x, abs(y) )
POW_EXP    x = pow( x, y )
POW_BASE   x = pow( y, x )
=========  =================================

The enum class is used in order to have one common function with an addtional op argument instead
of individual functions for each kind of operator.

