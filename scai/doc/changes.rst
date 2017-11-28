.. _changes:

Changes Version 3.0
===================

Here is a summary of the major changes:

 * The build system has been restructured and simplified. There are no more separate build steps 
   for each individual LAMA module. As a build of a module required the dependent modules to be
   installed before, build and installation was mixed up and lead to many confusions and inconveniences.
   Experience has also shown that individual settings for each module caused more troubles than 
   giving higher flexibility.
 * LAMA now uses full functionality of C++11, therefore the C++ compiler used for LAMA must
   support C++11 features.
 * The classes Matrix and Vector are now template classes with the value type of the elements 
   as template argument. There are still the classes _Matrix and _Vector as common base
   classes, but operations on these base classes are rather limited.
 * In contrary to the previous version, writing code with matrix-vector operations for arbitrary
   types is no longer supported. As a consequence, all solver classes are now also template classes.
 * The class Scalar as common container type for all kind of value types is no longer supported.
   Especially all methods on vector and scalars that return a single value return now a value of 
   the corresponding ValueType and no more a scalar. The class Scalar is still available as it is
   used internally for representation of scalar values in template expressions used to support text book
   syntax, but should not be used any more in user applications itself.
 * Enum types are now all enum classes (C++11 feature). Due to the scoping rules, the qualification of enum values
   is now slightly differrent.
 * The support for matrix-vector expressions has slightly improved and is now more functional.

Typed Vector Classes
--------------------

Here is a typical example of a flexible old LAMA code.

.. code-block:: c++

    void subroutine( Vector& x, const Matrix& m, const Vector& y )
    {
        VectorPtr z( y.newVector() );
        *z = m * y - x;
        Scalar alpha = z->l2Norm();
        x += alpha * y;
    }

With the latest version, this code has to be rewritten as follows:

.. code-block:: c++

    template<typename ValueType>
    void subroutine( Vector<ValueType>& x, const Matrix<ValueType>& m, const Vector<ValueType>& y )
    {
        DenseVector<ValueType> z;
        z = m * y - x;
        ValueType alpha = z.l2Norm();
        x += alpha * z;
    }

The is no more need to use a pointer variable for the temporary vector as now the type
of the vector is explicitly given. Furthermore the use of the class Scalar is no
more necessary and the used ValueType more explicit.

Matrix-Vector Expression
------------------------

.. code-block:: c++

     DenseVector<float> x( 100, 1 );
     DenseVector<double> y( 100, 2 );
     DenseVector<float> z;
     z = x + y;       // you cannot mix different types in expressions

Enum Types
----------

Example

.. code-block:: c++

    common::BinaryOp op = common::BinaryOp::COPY;  // common::binary::BinaryOp op = binary::COPY;

Matrix-Vector Expressions
-------------------------

.. code-block:: c++

   DenseVector<double> x;
   x.setRandom( 1000, 1 );    // 1000 random value between 0 and 1
   x = sin ( x );             // was x.sin() before
   x = cos ( x );             // was x.cos() before
   y = exp ( x );             // same as y = x; y.exp() before



