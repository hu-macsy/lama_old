.. _BinaryOp:

BinaryOp
========

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
POW        x = pow( x, y )
COPY_SIGN  x = copysign( x, y )
=========  =================================

The enum class is used in order to have one common function with an addtional op argument instead
of individual functions for each kind of operator.

The ``COPY`` operation is especially used in setter methods. Instead of combining old and
new values, the corresponding values are just replaced.

UnaryOp
=======

The following enumeration type specifies the different kind of unary operators
that can be used in elemenwise operations on arrays.

=========  =================================
Name       Operation
=========  =================================
CONJ       for conjugate 
ABS        for absolute value
MINUS      for negative value
EXP        call exp on each vector element
SQRT       call sqrt on each vector element
SIN        call sin on each vector element
COS        trigonometric function cos for each vector element
TAN        trigonometric function tan on each vector element
ATAN       call atan on each vector element
LOG        call log on each vector element
FLOOR      rounds downward
CEIL       rounds upward
=========  =================================

The enum class is used in order to have one common function with an addtional op argument instead
of individual functions for each kind of operator.
