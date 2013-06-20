Expressions
===========

Text Book Syntax in LAMA
------------------------

LAMA supports the use of text book syntax for expressions on matrices and vectors. Therefore, it is possible to use binary
operators like ``+``, ``-``, and ``*`` in your program when using LAMA matrices and vectors.

.. code-block:: c++

   #include <lama.hpp>
   #include <lama/Vector.hpp>
   #include <lama/matrix/all.hpp>

   void sub( Matrix& m, const Matrix& m1, const Matrix& m2, const Matrix& m3,
             Vector& v, const Vector& v1, const Vector& v2)
   {
       ...
       v = 1.5 * v2 + v1 * 2.5;    // vector expression
       v = 0.5 * m1 * v2 + s * v3; // matrix-vector expression
       m = 0.5 * m1 + m2;          // matrix expression
       m = 0.5 * m1 * m2 + m3;     // matrix-matrix expression
   }

The idea of matrix and vector expressions in LAMA is that these expressions are built in a first step as 'symbolic'
expressions and the symbolic expressions will not be evaluated before an assignment operator is called.
The assignment operator will select the corresponding factors of the 'symbolic' expression and call corresponding
methods to evaluate the full assignment.

.. code-block:: c++

   m = 0.5 * m1 * m2 + m3;  // ->  m1.matrixTimesMatrix( m, 0.5, m2, 1.0, m3)
   m = 0.5 * m1 + m2;       // ->  m.matrixPlusMatrix( 0.5, m1, 1.0, m2 )

By this strategy, LAMA avoids the use of intermediate vectors or matrices as results of binary operations; wherever it
is possible the whole expression is evaluated in one single call. Another advantage of these symbolic expressions are
that they can be normalized and optimized before they are finally resolved to corresponding operations.

.. code-block:: c++

   v = alpha * v1 * beta + gamma * v2; // v.vectorPlusVector( alpha * beta, v1, gamma, v2)

Unsupported expressions will already fail at compilation time.

.. code-block:: c++

   v = 0.5 * m1 * v2;  // OKAY, legal matrix-vector expression
   v = 0.5 * v2 * m1;  // ERROR, illegal matrix-vector expression
   v = v1 + v2 + v3;   // ERROR, unsupported as too complex

Size, Type, and Distribution Conformance
----------------------------------------

When using matrix and vector expressions, the operands must be conform, i.e. corresponding sizes of rows and columns 
must fit with other ones. LAMA will check these conditions, but this can only be done at runtime.

.. code-block:: c++

    CSRSparseMatrix<double> A( 3, 5), B( 5, 4 );
    CSRSparseMatrix<double> C = A * B;   // okay
    CSRSparseMatrix<double> D = B * A;   // runtime error

Multiprecision support has been considered as an important design goal of LAMA but is not yet implemented entirely.
Therefore, the operands should have all the same value type. In some situations LAMA supports implicit type conversions,
especially in simple assignments. But using different precisions in expressions is not always supported and will give
runtime errors.

.. code-block:: c++

    DenseVector<double> v1( 10, 5.0 );
    DenseVector<double> v2( 10, 5.0 );
    DenseVector<float>  v;

    v = v1 + v2;   // not supported due to type conversion, RUNTIME ERROR

    // Working solution

    DenseVector<double> tmp( v1 + v2 );
    v = tmp;       

Furthermore, the operands in the expressions should have the same distribution.

   - Distribution of vectors v1 and v2 must be the same in ``alpha * v1 + beta * v2``
   - Row and column distribution of matrices m1 and m2 must be the same in ``alpha * m1 + beta * m2``
   - Column distribution of m should be same as distribution of v in ``alpha * m * v``
   - Column distribution of m1 should be same as row distribution of m2 in ``alpha * m1 * m2``

The results of matrix and vector expressions will inherit distribution as the operands.

.. code-block:: c++

   // should/must be valid: v1.getDistribution() == v2.getDistribution()

   v = 1.5 * v2 + v1 * 2.5;    

   // now it is valid: v.getDistribution() == v1.getDistribution()

Unfortunately, it is not always easy to identify which expression has failed in the conformance checks. At least the
debug version of LAMA will print a call stack that might be very helpul to identify the source code line where the wrong
expression appears.

Vector-Expressions
------------------

A vector expression is a sum of scaled vectors. One summand can be used for incrementation or decrementation of a
vector, up to two summands are supported in an assignment.

.. code-block:: c++

   void sub( Vector& v, const Vector& v1, const Vector& v2, const Vector& v3 )
   Scalar alpha( 1.5 );
   Scalar beta( 2.0 );

   v = alpha * v1;
   v = alpha * v1 + beta * v2;
   v += alpha * v1;
   v -= alpha * v1;
 
The scalars used as scaling factors for the vectors will be represented as Scalar objects. Implicit type conversions
from double, int, or float values to Scalar are supported, so values of these types can be used in vector expressions
at any time.

When building symbolic vector expressions (Expression_SV), the following normalizations are done:

   * v1 becomes 1.0 * v1
   * v1 * alpha becomes alpha * v1
   * v1 / alpha becomes ( 1.0 / alpha )  * v1


Matrix-Vector-Expressions
-------------------------

A matrix-vector expression is a scaled matrix-vector product.

.. code-block:: c++

    void sub( Vector& v, const Matrix& m, const Vector v1 )
    Scalar alpha;
    ...
    v = alpha * m * v1;
    v += alpha * m * v1;
    v -= alpha * m * v1;

The size of vector ``v1`` must be equal to the number of columns in the matrix. The size of the result vector will be
equal to the number of rows of the matrix.

When building symbolic matrix-vector expressions (``Expression_SMV``) , the following normalizations are done:

   * ``m * v`` becomes ``1.0 * m * v``
   * ``m * v * alpha`` becomes ``alpha * m * v``
   * ``alpha * m * v * beta``  becomes ``(alpha * beta) * m  * v``
   * ``m * v / alpha``  becomes ``( 1.0 / alpha ) * m  * v``

A matrix-vector expression (``Expression_SMV``) can be added with a vector expression (``Expression_SV``)
and gives an expression (``Expression_SMV_SV``) also supported in an assignment.

.. code-block:: c++

    void sub( Vector& v, const Matrix& m, const Vector v1, const Vector v2 )
    {   
        Scalar alpha, beta
        ...
        v = alpha * m * v1 + beta * v2;

Matrix-Expressions
------------------

A matrix expression is a sum of scaled matrices. One summand can be used for incrementation or decrementation of
a matrix, up to two summands are supported in an assignment.

.. code-block:: c++

   void sub( Matrix& m, const Matrix& m1, const Matrix& m2 )
   Scalar alpha
   Scalar beta

   m = alpha * m1
   m = alpha * m1 + beta * m2
   m += alpha * m1
   m -= alpha * m1
 
When building symbolic matrx expressions, the following normalizations are done:

   * m1 becomes 1.0 * m1
   * m1 * alpha becomes alpha * m1
   * m1 / alpha becomes ( 1.0 / alpha )  * m1

In this sense, matrix expressions have nearly the same support as vector expressions.

Matrix-Matrix-Expression
------------------------

A matrix-matrix expression is a scaled matrix-matrix product.

.. code-block:: c++

    void sub( Matrix& m, const Matrix& m1, const Matrix& m2 )
    Scalar alpha;

    m = alpha * m1 * m2 ;
    m += alpha * m1 * m2;
    m -= alpha * m1 * m2;

When building symbolic matrix-matrix expressions, the following normalizations are done:

   * ``m1 * m2`` becomes ``1.0 * m1 * m2``
   * ``m1 * m2 * alpha`` becomes ``alpha * m1 * m2``
   * ``m1 * alpha * m2`` alpha becomes ``alpha * m1 * m2``

For the matrix-matrix product, the number of columns of the first matrix must be equal to the number of rows of the
second matrix. In case of distributed matrices, the column distribution of the first matrix should be equal to the row
distribution of the second matrix. It might be possible that LAMA can handle different distributions, but will at least
redistribute one of the matrices that might cause a certain overhead. For the result matrix, its row distribution will
be that of the first matrix, and its column distribution that of the second matrix.

In an assignment, a matrix-matrix expression can be added with a matrix expression.

.. code-block:: c++

    void sub( Matrix& m, const Matrix& m1, const Matrix& m2, const Matrix& m3 )
    Scalar alpha, beta;
    ...
    m = alpha * m1 * m2  + beta * m3;

Supported Expressions in Assignments
------------------------------------

.. code-block:: c++

    void sub( Vector& v, const Matrix& m, Vector& v1, const Vector& v2 )
    Scalar alpha, beta
    ...
    v = alpha * v1
    v = alpha * v1 + beta * v2
    v = alpha * m * v1 + beta * v2

The following expressions are supported for an assignment to a matrix:

   * matrix-expression ``scalar *  matrix``
   * matrix-expression ``scalar *  matrix + scalar * matrix``
   * matrix-expression ``scalar *  matrix * matrix``
   * matrix-expression ``scalar *  matrix * matrix + scalar * matrix``
   * and all expressions that can be transformed in such expressions

.. code-block:: c++

    void sub( Matrix& m, const Matrix& m1, const Matrix& m2, const Matrix& m3 )

    m = alpha * m1;
    m = alpha * m1 + beta * m2;
    m = alpha * m1 * m2;
    m = alpha * m1 * m2 + beta * m3;

    // other expressions are also supported if they can be normalized to above forms

    m = m1
    m = m1 + m2
    m = m1 * m2 * alpha + m3 
    m = m1 * alpha m2 * alpha + m3 

Other Expressions
=================

The operator ``*`` can be used to form the dotproduct of two vectors and will give
as a result a Scalar. The result as a Scalar might be used in other expressions
as well.

.. code-block:: c++

    DenseVector<double> x( 5, 1.0 );
    DenseVector<double> y( 5, 2.0 );
    Scalar s = x * y;
    DenseVector<double> z1 ( 3, 1.0 );
    DenseVector<double> z = x * y * z1;
    DenseVector<double> z = z1 * x * y; // RUNTIME error, z1, x are not conform

For the computation of the norm the Vector class provides corresponding methods.

.. code-block:: c++

    ...
    alpha = v1.l1Norm();   // l1 Norm
    alpha = v1.l2Norm();   // l2 Norm of a vector
    alpha = v1.maxNorm();  // max Norm of a vector

Another solution is the use of the norm classes.

.. code-block:: c++

    #include <lama/norm/all.hpp>
    ...
    L1Norm l1norm;
    L2Norm l2norm;
    MaxNorm maxnorm;
    alpha = l1norm( v );   // l1 Norm
    alpha = l2norm( v );   // l2 Norm of a vector
    alpha = maxnorm( v );  // max Norm of a vector

This solution is especially recommended when using different norms.

.. code-block:: c++

    void sub ( const Norm& norm )
    ...
    alpha = norm( v );   // calculate norm as required

    sub( L1Norm() );
    sub( L2Norm() );
    sub( MaxNorm() );

In future versions of LAMA, these norm classes are expected to deal with more general expressions that might avoid the
use of temporary vectors in case of differences.

.. code-block:: c++

    #include <lama/norm/all.hpp>
    ...
    const Norm& norm = ....
    alpha = norm( v1 - v2 ); // compute norm for vector difference
    alpha = norm( m1 - m2 ); // compute norm of matrix difference, elementwise

Constructors With Expressions
-----------------------------

All expressions that are supported in an assignment, can also be used in a constructor of a matrix.

.. code-block:: c++

    Matrix& m = ..., m1 = ..., m2 = ..., m3 = ...
    Vector& v = ..., v1 = ..., v2 = ..., v3 = ...
    Scalar alpha, beta

    DenseVector<ValueType> v( alpha * v1 + beta * v2 )
    DenseVector<ValueType> v( alpha * m1 * v1 + beta * v2 )
    CSRSparseMatrix<ValueType> m( alpha * m1 + beta * m2 )
    ELLSparseMatrix<ValueType> m( alpha * m1 * m2 + beta * m3 )

Performance Issues
------------------

Due to the use of symbolic expressions implememented by expression templates there is no performance loss for the
supported matrix and vector expressions. The little overhead is rather small and might be neglected for larger vectors
and matrices.

Regarding matrix and vector operations it is recommended that the operands have the same distribution. Even if LAMA
takes sometimes care of implicit redistributions, the corresponding overhead might slow down the performance.

