/**
 * @file ExpressionTest.cpp
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class ExpressionTest
 * @author Alexander Büchel, mdrost
 * @date 02.02.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>

#include <lama/CommunicationPlan.hpp>

#include <lama/Scalar.hpp>
#include <lama/DenseVector.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

struct ExpressionTestConfig
{
    ExpressionTestConfig()
    {
        numRows = 2;
        numCols = 2;
        rowDist = DistributionPtr( new NoDistribution( numRows ) );
        Dense.redistribute( rowDist, rowDist );
        Dense2.redistribute( rowDist, rowDist );
        Dense3.redistribute( rowDist, rowDist );
        x.resize( rowDist );
        x = 1.0;
        y.resize( rowDist );
        y = 0.0;
        alpha = 2.5;
        beta = -2.5;
    }

    ~ExpressionTestConfig()
    {
        rowDist = DistributionPtr();
    }

    IndexType numRows;
    IndexType numCols;
    DistributionPtr rowDist;
    DenseMatrix<double> Dense;
    DenseMatrix<double> Dense2;
    DenseMatrix<double> Dense3;
    DenseVector<double> x;
    DenseVector<double> y;
    Scalar alpha;
    Scalar beta;
};

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( ExpressionTest, ExpressionTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ExpressionTest" );

typedef Expression<Expression<Scalar,Expression<Matrix,Vector,Times>,Times>,Expression<Scalar,Vector,Times>,Plus> gemvExp;

typedef Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> gemmExp;

/* --------------------------------------------------------------------- */

void comparePointer( const void* p1, const void* p2, const int line, const std::string variable )
{
    std::string message( "failure at line: " );
    message += boost::lexical_cast<std::string,int>( line );
    message += " on variable: " + variable;
    BOOST_CHECK_MESSAGE( p1 == p2, message );
}

void compareScalar( const Scalar p1, const Scalar p2, const int line, const std::string variable )
{
    std::string message( "failure at line: " );
    message += boost::lexical_cast<std::string,int>( line );
    message += " on variable: " + variable;
    BOOST_CHECK_MESSAGE( p1 == p2, message );
}

void testVectorTimesScalar(
    const Scalar& alpha,
    const Vector& x,
    const Expression<Scalar,Vector,Times>& exp,
    const int line )
{
    compareScalar( exp.getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg2(), &x, line, "x" );
}

void testMatrixTimesScalar(
    const Scalar& alpha,
    const Matrix& m,
    const Expression<Scalar,Matrix,Times>& exp,
    const int line )
{
    compareScalar( exp.getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg2(), &m, line, "A" );
}

void testMatrixTimesScalarTimesVector(
    const Vector& x,
    const Matrix& m,
    const Scalar& alpha,
    const Expression<Scalar,Expression<Matrix,Vector,Times>,Times>& exp,
    const int line )
{
    compareScalar( exp.getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg2().getArg1(), &m, line, "A" );
    comparePointer( &exp.getArg2().getArg2(), &x, line, "x" );
}

void testMatrixTimesVector(
    const Matrix& m,
    const Vector& x,
    const Expression<Matrix,Vector,Times>& exp,
    const int line )
{
    comparePointer( &exp.getArg2(), &x, line, "x" );
    comparePointer( &exp.getArg1(), &m, line, "A" );
}

void testMatrixTimesMatrix(
    const Matrix& A,
    const Matrix& B,
    const Expression<Matrix,Matrix,Times>& exp,
    const int line )
{
    comparePointer( &exp.getArg1(), &A, line, "A" );
    comparePointer( &exp.getArg2(), &B, line, "B" );
}

void testScalarTimesMatrixTimesMatrix(
    const Scalar& alpha,
    const Matrix& A,
    const Matrix& B,
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp,
    const int line )
{
    compareScalar( exp.getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg2().getArg1(), &A, line, "A" );
    comparePointer( &exp.getArg2().getArg2(), &B, line, "B" );
}

void testGEMV(
    const Scalar& alpha,
    const Matrix& m,
    const Vector& x,
    const Scalar& beta,
    const Vector& y,
    const gemvExp& exp,
    const int line )
{
    compareScalar( exp.getArg1().getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg1().getArg2().getArg1(), &m, line, "A" );
    comparePointer( &exp.getArg1().getArg2().getArg2(), &x, line, "x" );
    compareScalar( exp.getArg2().getArg1(), beta, line, "beta" );
    comparePointer( &exp.getArg2().getArg2(), &y, line, "y" );
}

void testGEMM(
    const Scalar& alpha,
    const Matrix& A,
    const Matrix& B,
    const Scalar& beta,
    const Matrix& C,
    const gemmExp& exp,
    const int line )
{
    compareScalar( exp.getArg1().getArg1(), alpha, line, "alpha" );
    comparePointer( &exp.getArg1().getArg2().getArg1(), &A, line, "A" );
    comparePointer( &exp.getArg1().getArg2().getArg2(), &B, line, "B" );
    compareScalar( exp.getArg2().getArg1(), beta, line, "beta" );
    comparePointer( &exp.getArg2().getArg2(), &C, line, "C" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorScalar )
{
    testVectorTimesScalar( alpha, x, alpha * x, __LINE__ );

    testVectorTimesScalar( alpha, x, x * alpha, __LINE__ );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScalarTest )
{
    Scalar a( 1.0 );
    Scalar b( 2.0 );
    Scalar c;
    Scalar d( 4.0 );
    Scalar e( 5.0 );
    Scalar f;

    c = a + b;
    BOOST_CHECK( c == Scalar( 3.0 ) );

    c = b - a;
    BOOST_CHECK( c == Scalar( 1.0 ) );

    c = b * a;
    BOOST_CHECK( c == Scalar( 2.0 ) );

    c = a / b;
    BOOST_CHECK( c == Scalar( 0.5 ) );

    c = 4.0;
    c = sqrt( c );
    BOOST_CHECK( c == Scalar( 2.0 ) );

    c = 0.0;
    c = max( a, b );
    BOOST_CHECK( c == Scalar( 2.0 ) );

    c = -5.0;
    c = abs( c );
    BOOST_CHECK( c == Scalar( 5.0 ) );

    c = 3.0;
    f = ( b * c - a ) / e * d;
    BOOST_CHECK( f == Scalar( 4.0 ) );

    f = -100.0;
    f = ( b - a * c ) * e + d;
    BOOST_CHECK( f == Scalar( -1.0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MatrixScalar )
{
    testMatrixTimesScalar( alpha, Dense, alpha * Dense, __LINE__ );

    testMatrixTimesScalar( alpha, Dense, Dense * alpha, __LINE__ );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MatrixVector )
{
    testMatrixTimesVector( Dense, x, Dense * x, __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, alpha * ( Dense * x ), __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, ( Dense * x ) * alpha, __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, ( alpha * Dense ) * x, __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, Dense * ( x * alpha ), __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, alpha * Dense * x, __LINE__ );

    testMatrixTimesScalarTimesVector( x, Dense, alpha, Dense * x * alpha, __LINE__ );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( MatrixMatrix )
{
    testMatrixTimesMatrix( Dense, Dense2, Dense * Dense2, __LINE__ );

    testScalarTimesMatrixTimesMatrix( alpha, Dense, Dense2, alpha * ( Dense * Dense2 ), __LINE__ );

    testScalarTimesMatrixTimesMatrix( alpha, Dense, Dense2, ( alpha * Dense ) * Dense2, __LINE__ );

    testScalarTimesMatrixTimesMatrix( alpha, Dense, Dense2, Dense * ( alpha * Dense2 ), __LINE__ );

    testScalarTimesMatrixTimesMatrix( alpha, Dense, Dense2, Dense * alpha * Dense2, __LINE__ );

    testScalarTimesMatrixTimesMatrix( alpha, Dense, Dense2, alpha * Dense * Dense2, __LINE__ );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GEMV )
{
    testGEMV( alpha, Dense, x, beta, y, ( alpha * ( Dense * x ) ) + ( beta * y ), __LINE__ );

    testGEMV( alpha, Dense, x, beta, y, ( beta * y ) + ( alpha * ( Dense * x ) ), __LINE__ );

    // - --> beta*=-1 --> beta*-1=alpha
    testGEMV( alpha, Dense, x, alpha, y, ( alpha * ( Dense * x ) ) - ( beta * y ), __LINE__ );

    // - --> alpha*=-1 --> alpha*-1=beta
    testGEMV( beta, Dense, x, beta, y, ( beta * y ) - ( alpha * ( Dense * x ) ), __LINE__ );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GEMM )
{
    testGEMM( alpha, Dense, Dense2, beta, Dense3, alpha * Dense * Dense2 + beta * Dense3, __LINE__ );

    testGEMM( alpha, Dense, Dense2, beta, Dense3, beta * Dense3 + alpha * Dense * Dense2, __LINE__ );

    // - --> beta*=-1 --> beta*-1=alpha
    testGEMM( alpha, Dense, Dense2, alpha, Dense3, alpha * Dense * Dense2 - beta * Dense3, __LINE__ );

    // - --> alpha*=-1 --> alpha*-1=beta
    testGEMM( beta, Dense, Dense2, beta, Dense3, beta * Dense3 - alpha * Dense * Dense2, __LINE__ );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
