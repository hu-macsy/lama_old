/**
 * @file DenseMatrixTest.cpp
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
 * @brief Contains the implementation of the class DenseMatrixTest
 * @author Alexander Büchel
 * @date 03.04.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <test/TestSparseMatrices.hpp>
#include <test/TestMacros.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/expression/MatrixExpressions.hpp>

using namespace boost;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.DenseMatrixTest" );

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

template<typename T>
void verifyMatrixWithScalar( Matrix& m, Scalar& s )
{
    typedef T ValueType;

    for( IndexType i = 0; i < m.getNumRows(); ++i )
    {
        for( IndexType j = 0; j < m.getNumColumns(); ++j )
        {
            Scalar expected = m.getValue( i, j );
            Scalar result = s;
            LAMA_CHECK_SCALAR_SMALL_EPS( expected - result, ValueType );
        }
    }
}

template<typename T>
void verifySameMatrix( Matrix& m1, Matrix& m2 )
{
    typedef T ValueType;

    BOOST_REQUIRE_EQUAL( m1.getNumRows(), m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( m1.getNumColumns(), m2.getNumColumns() );

    for( IndexType i = 0; i < m1.getNumRows(); ++i )
    {
        for( IndexType j = 0; j < m1.getNumColumns(); ++j )
        {
            Scalar expected = m1.getValue( i, j );
            Scalar result = m2.getValue( i, j );
            LAMA_CHECK_SCALAR_SMALL_EPS( expected - result, ValueType );
        }
    }
}

template<typename T>
void verifySameMatrix( Matrix& m1, Matrix& m2, T eps )
{
    typedef T ValueType;

    BOOST_REQUIRE_EQUAL( m1.getNumRows(), m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( m1.getNumColumns(), m2.getNumColumns() );

    for( IndexType i = 0; i < m1.getNumRows(); ++i )
    {
        for( IndexType j = 0; j < m1.getNumColumns(); ++j )
        {
            Scalar expected = m1.getValue( i, j );
            Scalar result = m2.getValue( i, j );
            LAMA_LOG_INFO( logger,
                           "compare exp =" << expected.getValue<ValueType>() << ", res = " << result.getValue<ValueType>() );
            LAMA_CHECK_SCALAR_SMALL( expected - result, ValueType, eps );
        }
    }
}

/* -------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AssignmentMultiplicationTest, T, test_types ) {
    typedef T ValueType;
    IndexType n = 4;

//4x4 * 4x4
    {
        const Scalar s = 2.0;
        const Scalar t = 1.0;

        ValueType values[] =
        {   2.0f, 0.0f, 0.0f, 0.0f,
            0.0f, 2.0f, 0.0f, 0.0f,
            0.0f, 0.0f, 2.0f, 0.0f,
            0.0f, 0.0f, 0.0f, 2.0f
        };

        CSRSparseMatrix<ValueType> Diag2Matrix;
        Diag2Matrix.setRawDenseData( n, n, values );

        CSRSparseMatrix<ValueType> n4m4Matrix =
            TestSparseMatrices::n4m4TestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> n4m4InvMatrix =
            TestSparseMatrices::n4m4InverseTestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> identMatrix =
            TestSparseMatrices::n4m4IdentityMatrix<ValueType>();

        DenseMatrix<ValueType> matrixA( n4m4Matrix );
        const DenseMatrix<ValueType> matrixAInv( n4m4InvMatrix );
        DenseMatrix<ValueType> matrixIdent( identMatrix );

        DenseMatrix<ValueType> matrixD( 4, 4 );

        LAMA_LOG_INFO( logger, "linear algebra expression: a*A*B+b*C" );
        matrixD = t * matrixAInv * matrixA + t * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );

        LAMA_LOG_INFO( logger, "linear algebra expression: A*(A^-1)=Ident" );
        matrixD = matrixA * matrixAInv;
        verifySameMatrix<ValueType>( matrixD, matrixIdent );

        LAMA_LOG_INFO( logger, "selfassignment: A=A*B by hand: tmp=A, A=tmp*B" );
        DenseMatrix<ValueType> tmpA ( matrixA );
        matrixA = tmpA * matrixAInv;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );

        LAMA_LOG_INFO( logger, "selfassignment lhs - linear algebra expression: A=A*B" );
        matrixA = n4m4Matrix;
        matrixA = matrixA * matrixAInv;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );

        LAMA_LOG_INFO( logger, "reset and selfassignment rhs - linear algebra expression: A=B*A" );
        matrixA = n4m4Matrix;
        matrixA = matrixAInv * matrixA;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );

        LAMA_LOG_INFO( logger, "selfassignment lhs and rhs - linear algebra expression: A=A*A" );
        DenseMatrix<ValueType> matrixIdent2( identMatrix );
        matrixIdent = matrixIdent * matrixIdent;
        verifySameMatrix<ValueType>( matrixIdent2, matrixIdent );

        LAMA_LOG_INFO( logger, "linear algebra expression: a*A*B" );
        matrixA = n4m4Matrix;
        matrixD = s * matrixAInv * matrixA;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );

        LAMA_LOG_INFO( logger, "linear algebra expression: a*A" );
        matrixD = matrixA;
        matrixIdent = identMatrix;
        matrixD = s * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );

        LAMA_LOG_INFO( logger, "linear algebra expression: M*a" );
        matrixD = matrixA;
        matrixD = matrixIdent * s;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );

        LAMA_LOG_INFO( logger, "clear - linear algebra expression: a*A" );
        matrixD.clear();
        matrixD = s * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );

    } //4x4 * 4x4

//6x4 * 4x6, Constructor Test
    {
        CSRSparseMatrix<ValueType> n4m6SparseMatrix =
            TestSparseMatrices::n4m6MatrixD2<ValueType>();

        CSRSparseMatrix<ValueType> n6m4SparseMatrix =
            TestSparseMatrices::n6m4MatrixD1<ValueType>();

        CSRSparseMatrix<ValueType> resSparseMatrix =
            TestSparseMatrices::n6m6MatrixDRes<ValueType>();

        DenseMatrix<ValueType> n4m6DMatrix( n4m6SparseMatrix );
        DenseMatrix<ValueType> n6m4DMatrix( n6m4SparseMatrix );
        DenseMatrix<ValueType> resDMatrix( resSparseMatrix );

        DenseMatrix<ValueType> computeMatrix( n6m4DMatrix * n4m6DMatrix );

        verifySameMatrix<ValueType>( computeMatrix, resDMatrix );

        BOOST_CHECK( computeMatrix.getNumRows() == n6m4DMatrix.getNumRows() );
        BOOST_CHECK( computeMatrix.getNumColumns() == n4m6DMatrix.getNumColumns() );

        BOOST_CHECK( !(computeMatrix.getNumRows() == n4m6DMatrix.getNumRows() ) );
        BOOST_CHECK( !(computeMatrix.getNumColumns() == n6m4DMatrix.getNumColumns() ) );
    } //6x4 * 4x6

//ScalarMatrixMatrix
    {
        ValueType values[] =
        {   1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f
        };

        DenseMatrix<ValueType> mat1;
        mat1.setRawDenseData( 4, 4, values );
        DenseMatrix<ValueType> matrixRes( mat1 );
        ValueType j = 2.0;
        Scalar s = 2.0;
        Scalar t = 4.0;

        LAMA_LOG_INFO( logger, "matrixRes = j (" << j << ") * mat1 = " << mat1 );
        matrixRes = j * mat1;
        LAMA_LOG_INFO( logger, "matrixRes = " << matrixRes );
        verifyMatrixWithScalar<ValueType>( matrixRes, s );

        LAMA_LOG_INFO( logger, "matrixRes = mat1 ( " << mat1 << ") * j ( " << j << " )" );
        matrixRes = mat1 * j;
        LAMA_LOG_INFO( logger, "matrixRes = " << matrixRes );
        verifyMatrixWithScalar<ValueType>( matrixRes, s );

        LAMA_LOG_INFO( logger, "matrixRes = matrixRes ( " << matrixRes << ") * j ( " << j << " )" );
        matrixRes = matrixRes * j;
        LAMA_LOG_INFO( logger, "matrixRes = " << matrixRes );
        verifyMatrixWithScalar<ValueType>( matrixRes, t );
    }

    {
        CSRSparseMatrix<ValueType> n4m4Matrix =
            TestSparseMatrices::n4m4TestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> n4m4InvMatrix =
            TestSparseMatrices::n4m4InverseTestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> n4m6Matrix =
            TestSparseMatrices::n4m6MatrixD2<ValueType>();

        CSRSparseMatrix<ValueType> n6m4Matrix =
            TestSparseMatrices::n6m4MatrixD1<ValueType>();

        DenseMatrix<ValueType> matrixA( n4m4Matrix );
        DenseMatrix<ValueType> matrixInvA( n4m4InvMatrix );
        DenseMatrix<ValueType> matrixB = n4m6Matrix;
        DenseMatrix<ValueType> matrixC = n6m4Matrix;

        LAMA_LOG_INFO( logger, "Assignment Test - using constructor" );
        verifySameMatrix<ValueType>( n4m4Matrix, matrixA );

        LAMA_LOG_INFO( logger, "Assignment Test - using operator=" );
        verifySameMatrix<ValueType>( matrixB, n4m6Matrix );

        LAMA_LOG_INFO( logger, "Assignment Test - different matrix-sizes" );
        matrixA = matrixC;
        verifySameMatrix<ValueType>( matrixA, matrixC );
    }

}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( XGEMMOperationTest, T, test_types) {
    typedef T ValueType;

//alpha * A * B + beta * C
    {
        CSRSparseMatrix<ValueType> n4m4Matrix =
            TestSparseMatrices::n4m4TestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> n4m4MatrixInv =
            TestSparseMatrices::n4m4InverseTestMatrix1<ValueType>();

        CSRSparseMatrix<ValueType> identMatrix =
            TestSparseMatrices::n4m4IdentityMatrix<ValueType>();

        DenseMatrix<ValueType> matrixA( n4m4Matrix );
        DenseMatrix<ValueType> matrixAInv( n4m4MatrixInv );
        DenseMatrix<ValueType> matrixIdent( identMatrix );
        DenseMatrix<ValueType> matrixRes( identMatrix );
        DenseMatrix<ValueType> matrixTestRes( identMatrix );
        Scalar j = 2.0;

        matrixTestRes = j * j * matrixTestRes;

        LAMA_LOG_INFO( logger, "matrixRes = j * matrixA * matrixAInv + j * matrixIdent" );
        matrixRes = j * matrixA * matrixAInv + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-4f );

        matrixRes = matrixA;
        LAMA_LOG_INFO( logger, "matrixRes = j * matrixMatrixRes * matrixAInv + j * matrixIdent" );
        matrixRes = j * matrixRes * matrixAInv + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f );

        matrixRes = matrixAInv;
        LAMA_LOG_INFO( logger, "matrixRes = j * matrixA * matrixRes + j * matrixIdent" );
        matrixRes = j * matrixA * matrixRes + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f );

        matrixRes = matrixIdent;
        LAMA_LOG_INFO( logger, "matrixRes = j * matrixA * matrixAInv + j * matrixRes" );
        matrixRes = j * matrixA * matrixAInv + j * matrixRes;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f );

        matrixRes = matrixIdent;
        LAMA_LOG_INFO( logger, "matrixRes = j * matrixRes * matrixRes + j * matrixRes" );
        matrixRes = j * matrixRes * matrixRes + j * matrixRes;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f );
    }

    {
        CSRSparseMatrix<ValueType> n6m4SparseMatrix =
            TestSparseMatrices::n6m4MatrixE1<ValueType>();
        CSRSparseMatrix<ValueType> n4m3SparseMatrix =
            TestSparseMatrices::n4m3MatrixE2<ValueType>();
        CSRSparseMatrix<ValueType> resSparseMatrix =
            TestSparseMatrices::n6m3MatrixERes<ValueType>();

        DenseMatrix<ValueType> n6m4DMatrix( n6m4SparseMatrix );
        DenseMatrix<ValueType> n4m3DMatrix( n4m3SparseMatrix );
        DenseMatrix<ValueType> resDMatrix( 2.0 * resSparseMatrix );
        DenseMatrix<ValueType> cDMatrix( resSparseMatrix );

        DenseMatrix<ValueType> ergDMatrix( 6, 3 );

        Scalar j = 1.0;

        LAMA_LOG_INFO( logger, "ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix" );
        ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix;
        verifySameMatrix<ValueType>( ergDMatrix, resDMatrix );

        LAMA_LOG_INFO( logger, "cDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix" );
        cDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix;
        verifySameMatrix<ValueType>( ergDMatrix, resDMatrix );
    } //alpha * A * B + beta * C

//alpha * A * B
    {
        CSRSparseMatrix<ValueType> n4m4SMatrixA =
            TestSparseMatrices::n4m4MatrixA1<ValueType>();

        CSRSparseMatrix<ValueType> n4m4SMatrixB =
            TestSparseMatrices::n4m4MatrixA2<ValueType>();

        CSRSparseMatrix<ValueType> resSMatrix =
            TestSparseMatrices::n4m4MatrixARes<ValueType>();

        DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
        DenseMatrix<ValueType> matrixB( n4m4SMatrixB );
        DenseMatrix<ValueType> matrixRes( 4, 4 );
        DenseMatrix<ValueType> matrixTestRes( resSMatrix );
        ValueType j = 2.0;
        matrixTestRes = j * matrixTestRes;

        LAMA_LOG_INFO( logger, "matrixRes = j * matrixA * matrixB" );
        matrixRes = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes );

        LAMA_LOG_INFO( logger, "matrixA = j * matrixA * matrixB" );
        matrixA = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixA );

        matrixA = n4m4SMatrixA;
        LAMA_LOG_INFO( logger, "matrixB = j * matrixA * matrixB" );
        matrixB = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixB );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CtorTest, T, test_types ) {
    typedef T ValueType;

// C = A * B
//6x4 * 4x6, Constructor Test
    {
        CSRSparseMatrix<ValueType> n4m6SparseMatrix =
            TestSparseMatrices::n4m6MatrixD2<ValueType>();
        CSRSparseMatrix<ValueType> n6m4SparseMatrix =
            TestSparseMatrices::n6m4MatrixD1<ValueType>();
        CSRSparseMatrix<ValueType> resSparseMatrix =
            TestSparseMatrices::n6m6MatrixDRes<ValueType>();
        CSRSparseMatrix<ValueType> cSparseMatrix =
            TestSparseMatrices::n6m6TestMatrix<ValueType>();

        DenseMatrix<ValueType> n4m6DMatrix( n4m6SparseMatrix );
        DenseMatrix<ValueType> n6m4DMatrix( n6m4SparseMatrix );
        DenseMatrix<ValueType> resDMatrix( resSparseMatrix );

        BOOST_MESSAGE("TODO: DenseMatrixTest - CtorTest");
        DenseMatrix<ValueType> computeMatrix( n6m4DMatrix * n4m6DMatrix );

        verifySameMatrix<ValueType>( computeMatrix, resDMatrix );

        BOOST_CHECK( computeMatrix.getNumRows() == n6m4DMatrix.getNumRows() );
        BOOST_CHECK( computeMatrix.getNumColumns() == n4m6DMatrix.getNumColumns() );

        BOOST_CHECK( !(computeMatrix.getNumRows() == n4m6DMatrix.getNumRows() ) );
        BOOST_CHECK( !(computeMatrix.getNumColumns() == n6m4DMatrix.getNumColumns() ) );
    } //6x4 * 4x6

//a*A & A*a
    {
        Scalar j = 2.0;

        ValueType values[] =
        {   1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f,
            1.0f, 1.0f, 1.0f, 1.0f
        };

        DenseMatrix<ValueType> mat1;
        mat1.setRawDenseData( 4, 4, values );
        DenseMatrix<ValueType> matrixRes( j * mat1 );
        verifyMatrixWithScalar<ValueType>( matrixRes, j );

        DenseMatrix<ValueType> matrixRes2( mat1 * j );
        verifyMatrixWithScalar<ValueType>( matrixRes2, j );

    } //a*A & A*a

//    alpha*A*B+beta*C
//    XGEMM
    {
        CSRSparseMatrix<ValueType> n6m4SparseMatrix =
            TestSparseMatrices::n6m4MatrixE1<ValueType>();
        CSRSparseMatrix<ValueType> n4m3SparseMatrix =
            TestSparseMatrices::n4m3MatrixE2<ValueType>();
        CSRSparseMatrix<ValueType> resSparseMatrix =
            TestSparseMatrices::n6m3MatrixERes<ValueType>();

        DenseMatrix<ValueType>
        n6m4DMatrix( n6m4SparseMatrix );
        DenseMatrix<ValueType>
        n4m3DMatrix( n4m3SparseMatrix );
        DenseMatrix<ValueType>
        resDMatrix( 2.0 * resSparseMatrix );
        DenseMatrix<ValueType>
        cDMatrix( resSparseMatrix );

        Scalar j = 1.0;

        DenseMatrix<ValueType> ergDMatrix( 6,3 );
        ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix;
        verifySameMatrix<ValueType>( ergDMatrix, resDMatrix );

    } //alpha * A * B + beta * C

//alpha * A * B
    {
        CSRSparseMatrix<ValueType> n4m4SMatrixA =
            TestSparseMatrices::n4m4MatrixA1<ValueType>();

        CSRSparseMatrix<ValueType> n4m4SMatrixB =
            TestSparseMatrices::n4m4MatrixA2<ValueType>();

        CSRSparseMatrix<ValueType> resSMatrix =
            TestSparseMatrices::n4m4MatrixARes<ValueType>();

        DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
        DenseMatrix<ValueType> matrixB( n4m4SMatrixB );
        DenseMatrix<ValueType> matrixTestRes( resSMatrix );
        ValueType j = 2.0;
        matrixTestRes = j * matrixTestRes;

        DenseMatrix<ValueType> matrixRes( j * matrixA * matrixB );
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes );
    } //alpha * A * B
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetIdentityTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    matrixA.setIdentity();

    for( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        for( IndexType j = 0; j < matrixA.getNumColumns(); ++j )
        {
            if( i == j )
            {
                BOOST_CHECK_EQUAL( matrixA.getValue( i, j ), 1.0 );
            }
            else
            {
                BOOST_CHECK_EQUAL( matrixA.getValue( i, j ), 0.0 );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetDiagonalTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    DenseVector<ValueType> vector( 4, -1.0 );

    matrixA.setDiagonal( vector );

    for( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        BOOST_CHECK_EQUAL( matrixA.getValue( i, i ), -1.0 );
    }

    Scalar s = 4.0;

    matrixA.setDiagonal( s );

    for( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        BOOST_CHECK_EQUAL( matrixA.getValue( i, i ), 4.0 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetRowTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    DenseVector<ValueType> dvector( 0, 0 );

    Vector& vector = dvector;

    for( int i = 0; i < matrixA.getNumRows(); i++ )
    {
        matrixA.getRow( vector, i );

        for( int j = 0; j < matrixA.getNumColumns(); j++ )
        {
            BOOST_CHECK_EQUAL( vector.getValue( j ), matrixA.getValue( i, j ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetNumValuesTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    BOOST_CHECK_EQUAL( matrixA.getNumValues(), 9 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetUsageMemoryTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    size_t size = 141;

    BOOST_CHECK_EQUAL( matrixA.getMemoryUsage(), size );

}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetLocalValuesTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );

    BOOST_CHECK_EQUAL( matrixA.getLocalNumValues(), 16 );
    BOOST_CHECK_EQUAL( matrixA.getLocalNumRows(), 4 );
    BOOST_CHECK_EQUAL( matrixA.getLocalNumRows(), 4 );

}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScaleTest )
{
    typedef double ValueType;

    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();

    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    Scalar s = 2.0;
    matrixA.scale( s );

    for( int i = 0; i < matrixA.getNumRows(); i++ )
    {
        for( int j = 0; j < matrixA.getNumColumns(); j++ )
        {
            BOOST_CHECK_EQUAL( n4m4SMatrixA.getValue( i, j ) * 2.0, matrixA.getValue( i, j ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( TypeNameTest )
{
    CSRSparseMatrix<float> n4m4SMatrixAf = TestSparseMatrices::n4m4MatrixA1<float>();

    DenseMatrix<float> matrixA( n4m4SMatrixAf );
    const char* floatchar = matrixA.getTypeName();
    BOOST_CHECK_EQUAL( floatchar, "DenseMatrix<float>" );

    CSRSparseMatrix<double> n4m4SMatrixAd = TestSparseMatrices::n4m4MatrixA1<double>();

    DenseMatrix<double> matrixB( n4m4SMatrixAd );
    const char* doublechar = matrixB.getTypeName();
    BOOST_CHECK_EQUAL( doublechar, "DenseMatrix<double>" );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SwapTest, T, test_types ) {
    typedef T ValueType;

    CSRSparseMatrix<ValueType> SMatrixA =
        TestSparseMatrices::n4m4MatrixA1<ValueType>();

    CSRSparseMatrix<ValueType> SMatrixB =
        TestSparseMatrices::n4m4MatrixA2<ValueType>();

    DenseMatrix<ValueType> matrixA1( SMatrixA );
    DenseMatrix<ValueType> matrixB1( SMatrixB );

    DenseMatrix<ValueType> matrixA2( SMatrixA );
    DenseMatrix<ValueType> matrixB2( SMatrixB );

    matrixA1.swap( matrixB1 );

    verifySameMatrix<ValueType>( matrixA1, matrixB2 );
    verifySameMatrix<ValueType>( matrixB1, matrixA2 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
