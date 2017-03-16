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

.. code-block:: c++

    scai::common::binary::BinaryOp op = scai::common::binary::SUB;
    double x = 1.0;
    double y = 1.5;
    double res = scai::common::applyBinary<double>( x, op, y );

CompareOp
=========

The following enumeration type specifies the different kind of compare operators that
might be used.

=========  =================================
Name       Operation
=========  =================================
LE         x <= y
LT         x < y
GE         x >= y
GT         x > y
=========  =================================

.. code-block:: c++

    scai::common::binary::CompareOp op = scai::common::binary::LE;
    double x = 1.0;
    double y = 1.5;
    bool res = scai::common::applyBinary<double>( x, op, y );

UnaryOp
=======

The following enumeration type specifies the different kind of unary operators
that can be used in elemenwise operations on arrays.

=========  =================================
Name       Operation
=========  =================================
CONJ       for conjugate 
ABS        for absolute value
ASUM       same as ABS( real(x) ) + ABS( imag(x) )
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

.. code-block:: c++

    scai::common::unary::UnaryOp op = scai::common::unary::SIN;
    double x = 1.0;
    double y = scai::common::applyUnary<double>( op, x );

