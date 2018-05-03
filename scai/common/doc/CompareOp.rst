.. _CompareOp:

CompareOp
=========

The enumeration class ``CompareOp`` specifies the different kind of compare operators that
might be used. Its enumerators are:

=========  =================================
name       operation
=========  =================================
``LE``     x <= y
``LT``     x < y
``GE``     x >= y
``GT``     x > y
``EQ``     x == y
``NE``     x != y
=========  =================================

.. code-block:: c++

    using namespace scai;

    common::CompareOp op = ...
    std::cout << "comparison is " << op << std::endl;

    common::applyBinary<double>( 3.0, common::CompareOp::GT, 5.0 )  // returns false

    sort( ..., common::CompareOp::LE )          // sort ascending
    isSorted( ..., common:ComapreOp::GT );     // check for strong descedindg

    
