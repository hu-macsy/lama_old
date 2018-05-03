.. _UnaryOp:

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

    scai::common::UnaryOp op = scai::common::UnaryOp::SIN;
    double x = 1.0;
    double y = scai::common::applyUnary<double>( op, x );

