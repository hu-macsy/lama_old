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

    scai::common::BinaryOp op = scai::common::BinaryOp::SUB;
    double x = 1.0;
    double y = 1.5;
    double res = scai::common::applyBinary<double>( x, op, y );
