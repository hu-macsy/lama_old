/**
 * @file DenseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestSparseMatrices.hpp>
#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/SameMatrixHelper.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/all.hpp>

//-------------------------
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp> 
#include <scai/common/shared_ptr.hpp>   
//-------------------------
#include <scai/lama/matutils/MatrixCreator.hpp> 




using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using scai::common::unique_ptr;
using scai::common::scoped_array;
using scai::common::shared_ptr;
namespace scai
{
namespace lama
{
namespace DenseMatrixTest
{

/* ----------------------------- help functions --------------------------------------------------------------------- */

template<typename ValueType>
void verifyMatrixWithScalar( Matrix& m, Scalar& s )
{
    for ( IndexType i = 0; i < m.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < m.getNumColumns(); ++j )
        {
            Scalar expected = m.getValue( i, j );
            Scalar result = s;
            SCAI_CHECK_SCALAR_SMALL_EPS( expected - result, ValueType );
        }
    }
}

template<typename ValueType>
void verifySameMatrix( Matrix& m1, Matrix& m2 )
{
    BOOST_REQUIRE_EQUAL( m1.getNumRows(), m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( m1.getNumColumns(), m2.getNumColumns() );

    for ( IndexType i = 0; i < m1.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < m1.getNumColumns(); ++j )
        {
            Scalar expected = m1.getValue( i, j );
            Scalar result = m2.getValue( i, j );
            SCAI_CHECK_SCALAR_SMALL_EPS( expected - result, ValueType );
        }
    }
}

template<typename ValueType>
void verifySameMatrix( Matrix& m1, Matrix& m2, ValueType eps, logging::Logger& logger )
{
    BOOST_REQUIRE_EQUAL( m1.getNumRows(), m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( m1.getNumColumns(), m2.getNumColumns() );

    for ( IndexType i = 0; i < m1.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < m1.getNumColumns(); ++j )
        {
            Scalar expected = m1.getValue( i, j );
            Scalar result = m2.getValue( i, j );
            SCAI_LOG_INFO( logger,
                           "compare exp =" << expected.getValue<ValueType>() << ", res = " << result.getValue<ValueType>() );
            SCAI_CHECK_SCALAR_SMALL( expected - result, ValueType, eps );
        }
    }
}

template<typename ValueType>
void GEMMTestImpl( const int n, const int m, const int k, ValueType eps, ContextPtr loc )
{
    Scalar alpha( 2.0 );
    Scalar beta( 3.0 );
    int maxdim = std::max( n, m );
    // mainly for the very big example matrices the single values have to be chosen small
    // for the result not becoming too big for float!!!
    scoped_array<ValueType> values( new ValueType[maxdim * k] );

    for ( int i = 0; i < maxdim * k; i++ )
    {
        values[i] = i * static_cast<ValueType>( 1e-5 ) + 1;
    }

    scoped_array<ValueType> valuesC( new ValueType[m * n] );
    ValueType help;

    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            help = 0.0;

            for ( int kk = 0; kk < k; kk++ )
            {
                // stores row i (of data same like matrix A) times column k (of data same like matrix B)
                help += values[i * k + kk] * values[kk * n + j];
            }

            valuesC[i * n + j] = ( static_cast<ValueType>( 42.0 ) - alpha.getValue<ValueType>() * help ) / beta.getValue<ValueType>();
        }
    }

    DenseMatrix<ValueType> A;
    A.setContextPtr( loc, loc );
    A.setRawDenseData( m, k, values.get() );
    DenseMatrix<ValueType> B; // not possible: const B
    B.setRawDenseData( k, n, values.get() );
    DenseMatrix<ValueType> C;
    C.setRawDenseData( m, n, valuesC.get() );
    C = alpha * A * B + beta * C;

    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            Scalar expectedValue( 42.0 );
            Scalar value = C.getValue( i, j );
            Scalar diff = expectedValue - value;
            SCAI_CHECK_SCALAR_SMALL( diff, ValueType, eps );
        }
    }

    DenseMatrix<ValueType> C2;
    C2.setRawDenseData( m, n, valuesC.get() );
    DenseMatrix<ValueType> D( alpha * A * B + beta * C2 );

    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            Scalar expectedValue( 42.0 );
            Scalar value = D.getValue( i, j );
            Scalar diff = expectedValue - value;
            SCAI_CHECK_SCALAR_SMALL( diff, ValueType, eps );
        }
    }
}

/* ----------------------------- test functions --------------------------------------------------------------------- */

template<typename ValueType>
void assignmentMultiplicationTest( logging::Logger& logger )
{
    IndexType n = 4;
    //4x4 * 4x4
    {
        const Scalar s = 2.0;
        const Scalar t = 1.0;
        ValueType values[] =
        {
            2.0f, 0.0f, 0.0f, 0.0f,
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
        SCAI_LOG_INFO( logger, "linear algebra expression: a*A*B+b*C" );
        matrixD = t * matrixAInv * matrixA + t * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );
        SCAI_LOG_INFO( logger, "linear algebra expression: A*(A^-1)=Ident" );
        matrixD = matrixA * matrixAInv;
        verifySameMatrix<ValueType>( matrixD, matrixIdent );
        SCAI_LOG_INFO( logger, "selfassignment: A=A*B by hand: tmp=A, A=tmp*B" );
        DenseMatrix<ValueType> tmpA ( matrixA );
        matrixA = tmpA * matrixAInv;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );
        SCAI_LOG_INFO( logger, "selfassignment lhs - linear algebra expression: A=A*B" );
        matrixA = n4m4Matrix;
        matrixA = matrixA * matrixAInv;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );
        SCAI_LOG_INFO( logger, "reset and selfassignment rhs - linear algebra expression: A=B*A" );
        matrixA = n4m4Matrix;
        matrixA = matrixAInv * matrixA;
        verifySameMatrix<ValueType>( matrixA, matrixIdent );
        SCAI_LOG_INFO( logger, "self assignment lhs and rhs - linear algebra expression: A=A*A" );
        DenseMatrix<ValueType> matrixIdent2( identMatrix );
        matrixIdent = matrixIdent * matrixIdent;
        verifySameMatrix<ValueType>( matrixIdent2, matrixIdent );
        SCAI_LOG_INFO( logger, "linear algebra expression: a*A*B" );
        matrixA = n4m4Matrix;
        matrixD = s * matrixAInv * matrixA;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );
        SCAI_LOG_INFO( logger, "linear algebra expression: a*A" );
        matrixD = matrixA;
        matrixIdent = identMatrix;
        matrixD = s * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );
        SCAI_LOG_INFO( logger, "linear algebra expression: M*a" );
        matrixD = matrixA;
        matrixD = matrixIdent * s;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );
        SCAI_LOG_INFO( logger, "clear - linear algebra expression: a*A" );
        matrixD.clear();
        matrixD = s * matrixIdent;
        verifySameMatrix<ValueType>( matrixD, Diag2Matrix );
    } // 4x4 * 4x4
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
        BOOST_CHECK( !( computeMatrix.getNumRows() == n4m6DMatrix.getNumRows() ) );
        BOOST_CHECK( !( computeMatrix.getNumColumns() == n6m4DMatrix.getNumColumns() ) );
    } //6x4 * 4x6
    //ScalarMatrixMatrix
    {
        ValueType values[] =
        {
            1.0f, 1.0f, 1.0f, 1.0f,
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
        SCAI_LOG_INFO( logger, "matrixRes = j (" << j << ") * mat1 = " << mat1 );
        matrixRes = j * mat1;
        SCAI_LOG_INFO( logger, "matrixRes = " << matrixRes );
        verifyMatrixWithScalar<ValueType>( matrixRes, s );
        SCAI_LOG_INFO( logger, "matrixRes = mat1 ( " << mat1 << ") * j ( " << j << " )" );
        matrixRes = mat1 * j;
        SCAI_LOG_INFO( logger, "matrixRes = " << matrixRes );
        verifyMatrixWithScalar<ValueType>( matrixRes, s );
        SCAI_LOG_INFO( logger, "matrixRes = matrixRes ( " << matrixRes << ") * j ( " << j << " )" );
        matrixRes = matrixRes * j;
        SCAI_LOG_INFO( logger, "matrixRes = " << matrixRes );
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
        SCAI_LOG_INFO( logger, "Assignment Test - using constructor" );
        verifySameMatrix<ValueType>( n4m4Matrix, matrixA );
        SCAI_LOG_INFO( logger, "Assignment Test - using operator=" );
        verifySameMatrix<ValueType>( matrixB, n4m6Matrix );
        SCAI_LOG_INFO( logger, "Assignment Test - different matrix-sizes" );
        matrixA = matrixC;
        verifySameMatrix<ValueType>( matrixA, matrixC );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void xGEMMOperationTest( logging::Logger& logger )
{
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
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixA * matrixAInv + j * matrixIdent" );
        matrixRes = j * matrixA * matrixAInv + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-4f, logger );
        matrixRes = matrixA;
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixMatrixRes * matrixAInv + j * matrixIdent" );
        matrixRes = j * matrixRes * matrixAInv + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f, logger );
        matrixRes = matrixAInv;
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixA * matrixRes + j * matrixIdent" );
        matrixRes = j * matrixA * matrixRes + j * matrixIdent;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f, logger );
        matrixRes = matrixIdent;
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixA * matrixAInv + j * matrixRes" );
        matrixRes = j * matrixA * matrixAInv + j * matrixRes;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f, logger );
        matrixRes = matrixIdent;
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixRes * matrixRes + j * matrixRes" );
        matrixRes = j * matrixRes * matrixRes + j * matrixRes;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes, 1E-3f, logger );
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
        SCAI_LOG_INFO( logger, "ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix" );
        ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix;
        verifySameMatrix<ValueType>( ergDMatrix, resDMatrix );
        SCAI_LOG_INFO( logger, "cDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix" );
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
        SCAI_LOG_INFO( logger, "matrixRes = j * matrixA * matrixB" );
        matrixRes = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixRes );
        SCAI_LOG_INFO( logger, "matrixA = j * matrixA * matrixB" );
        matrixA = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixA );
        matrixA = n4m4SMatrixA;
        SCAI_LOG_INFO( logger, "matrixB = j * matrixA * matrixB" );
        matrixB = j * matrixA * matrixB;
        verifySameMatrix<ValueType>( matrixTestRes, matrixB );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void cTorTest( )
{
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
        BOOST_MESSAGE( "TODO: DenseMatrixTest - CtorTest" );
        DenseMatrix<ValueType> computeMatrix( n6m4DMatrix * n4m6DMatrix );
        verifySameMatrix<ValueType>( computeMatrix, resDMatrix );
        BOOST_CHECK( computeMatrix.getNumRows() == n6m4DMatrix.getNumRows() );
        BOOST_CHECK( computeMatrix.getNumColumns() == n4m6DMatrix.getNumColumns() );
        BOOST_CHECK( !( computeMatrix.getNumRows() == n4m6DMatrix.getNumRows() ) );
        BOOST_CHECK( !( computeMatrix.getNumColumns() == n6m4DMatrix.getNumColumns() ) );
    } //6x4 * 4x6
    //a*A & A*a
    {
        Scalar j = 2.0;
        ValueType values[] =
        {
            1.0f, 1.0f, 1.0f, 1.0f,
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
    } // a*A & A*a
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
        DenseMatrix<ValueType> ergDMatrix( 6, 3 );
        ergDMatrix = j * n6m4DMatrix * n4m3DMatrix + j * cDMatrix;
        verifySameMatrix<ValueType>( ergDMatrix, resDMatrix );
    } // alpha * A * B + beta * C
    // alpha * A * B
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
    } // alpha * A * B
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setIdentityTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA;
    matrixA.setIdentity( n4m4SMatrixA.getNumRows() );

    for ( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < matrixA.getNumColumns(); ++j )
        {
            if ( i == j )
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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setDiagonalTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    DenseVector<ValueType> vector( 4, -1.0 );
    matrixA.setDiagonal( vector );

    for ( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        BOOST_CHECK_EQUAL( matrixA.getValue( i, i ), -1.0 );
    }

    Scalar s = 4.0;
    matrixA.setDiagonal( s );

    for ( IndexType i = 0; i < matrixA.getNumRows(); ++i )
    {
        BOOST_CHECK_EQUAL( matrixA.getValue( i, i ), 4.0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void getRowTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    DenseVector<ValueType> dvector( 0, 0 );
    Vector& vector = dvector;

    for ( int i = 0; i < matrixA.getNumRows(); i++ )
    {
        matrixA.getRow( vector, i );

        for ( int j = 0; j < matrixA.getNumColumns(); j++ )
        {
            BOOST_CHECK_EQUAL( vector.getValue( j ), matrixA.getValue( i, j ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void getNumValuesTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    BOOST_CHECK_EQUAL( matrixA.getNumValues(), 9 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void getMemoryUsageTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    size_t size_float = 77;
    size_t size_double = 141;
    std::ostringstream omsg;
 
    omsg << common::getScalarType<ValueType>();

    if ( std::string( "double" ).compare( omsg.str() ) == 0 )
    {
        BOOST_CHECK_EQUAL( matrixA.getMemoryUsage(), size_double );
    }
    else
    {
        BOOST_CHECK_EQUAL( matrixA.getMemoryUsage(), size_float );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void getLocalValuesTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    BOOST_CHECK_EQUAL( matrixA.getLocalNumValues(), 16 );
    BOOST_CHECK_EQUAL( matrixA.getLocalNumRows(), 4 );
    BOOST_CHECK_EQUAL( matrixA.getLocalNumRows(), 4 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void scaleTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixA = TestSparseMatrices::n4m4MatrixA1<ValueType>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixA );
    Scalar s = 2.0;
    matrixA.scale( s );

    for ( int i = 0; i < matrixA.getNumRows(); i++ )
    {
        for ( int j = 0; j < matrixA.getNumColumns(); j++ )
        {
            BOOST_CHECK_EQUAL( n4m4SMatrixA.getValue( i, j ) * 2.0, matrixA.getValue( i, j ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void typeNameTest( )
{
    CSRSparseMatrix<ValueType> n4m4SMatrixAf = TestSparseMatrices::n4m4MatrixA1<float>();
    DenseMatrix<ValueType> matrixA( n4m4SMatrixAf );
    std::string s = matrixA.typeName();
    BOOST_CHECK( s.length() > 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void matrixAddMatrixTest( )
{
    DenseMatrix<ValueType> n6m4SMatrixA = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    DenseMatrix<ValueType> n6m4SMatrixB = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    // TODO ThoBra: Wenn man den Speicherplatz hier im Test nicht explizit alloziert, wird kein Speicherplatz
    // für mData bereitgestellt und memory access violation ist die Folge.
    DenseMatrix<ValueType> matrix( n6m4SMatrixA.getNumRows(), n6m4SMatrixA.getNumColumns() );
    matrix = n6m4SMatrixA + n6m4SMatrixB;
    DenseMatrix<ValueType> matrixResult = 2.0 * n6m4SMatrixA;
    testSameMatrix( matrix, matrixResult );
}

/* ------------------------------------------------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void maxNormTest( )
{
    DenseMatrix<ValueType> n6m4SMatrix = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    Scalar maxNorm = n6m4SMatrix.maxNorm();
    BOOST_CHECK_EQUAL( maxNorm.getValue<ValueType>(), 9.0 );
    DenseMatrix<ValueType> n6m6SMatrix = TestSparseMatrices::n6m6MatrixDRes<ValueType>();
    maxNorm = n6m6SMatrix.maxNorm();
    BOOST_CHECK_EQUAL( maxNorm.getValue<ValueType>(), 89.0 );
    DenseMatrix<ValueType> n6m4SMatrix1 = 2.5 * n6m4SMatrix;
    Scalar maxDiffNorm = n6m4SMatrix.maxDiffNorm( n6m4SMatrix1 );
    BOOST_CHECK_EQUAL( maxDiffNorm.getValue<ValueType>(), 13.5 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void gemmTest( ContextPtr loc )
{
    // This test needs less precision eps for float than for double
    ValueType eps;
    std::ostringstream omsg;
    omsg << common::getScalarType<ValueType>();

    if ( std::string( "double" ).compare( omsg.str() ) == 0 )
    {
        eps = static_cast<ValueType>( 1e-5 );
    }
    else
    {
        eps = static_cast<ValueType>( 1e-4 );
    }

    GEMMTestImpl<ValueType>( 5, 3, 4, eps, loc );
    GEMMTestImpl<ValueType>( 2, 2, 3, eps, loc );
    GEMMTestImpl<ValueType>( 16, 16, 16, eps, loc );
    GEMMTestImpl<ValueType>( 32, 16, 16, eps, loc );
    GEMMTestImpl<ValueType>( 2, 2, 4, eps, loc );
    GEMMTestImpl<ValueType>( 12, 12, 17, eps, loc );
    GEMMTestImpl<ValueType>( 32, 16, 32, eps, loc );
    GEMMTestImpl<ValueType>( 16, 32, 16, eps, loc );
    GEMMTestImpl<ValueType>( 32, 32, 16, eps, loc );
    GEMMTestImpl<ValueType>( 16, 32, 32, eps, loc );

    if ( std::string( "double" ).compare( omsg.str() ) == 0 )
    {
        eps = static_cast<ValueType>( 1e-4 );
    }
    else
    {
        eps = static_cast<ValueType>( 1e-2 );
    }

    GEMMTestImpl<ValueType>( 31, 31, 1114, eps, loc );
    GEMMTestImpl<ValueType>( 32, 32, 256, eps, loc );
    GEMMTestImpl<ValueType>( 64, 64, 64, eps, loc );
    GEMMTestImpl<ValueType>( 128, 128, 256, eps, loc );
}

/* ------------------------------------------------------------------------------------------------------------------ */
template<typename ValueType>
void swapTest( )
{
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



template<typename ValueType>
void transposeTest( )
{
    CommunicatorPtr comm = Communicator::getCommunicator(scai::dmemo::communicator::MPI);  // MPI
    CSRSparseMatrix<ValueType> testMatrixA = TestSparseMatrices::n6m4MatrixD1<ValueType>();
    CSRSparseMatrix<ValueType> testMatrixB = TestSparseMatrices::n6m6Full<ValueType>();
    DenseMatrix<ValueType> matrixA(testMatrixA);
    DenseMatrix<ValueType> matrixB(testMatrixB);
   
    // Block Distribution
    
    DistributionPtr rowsBlock( new BlockDistribution( matrixA.getNumRows(), comm ) );
    DistributionPtr columnsBlock(new BlockDistribution( matrixA.getNumColumns(), comm ) );
    matrixA.redistribute(rowsBlock,columnsBlock);
    
    DistributionPtr rowBlock( new BlockDistribution( matrixB.getNumRows(), comm ) );
    DistributionPtr columnBlock(new BlockDistribution( matrixB.getNumColumns(), comm ) );
    matrixB.redistribute(rowBlock,columnBlock);

    //Transpose
    matrixA.assignTranspose(matrixA);
    matrixB.assignTranspose(matrixB);
    matrixA.assignTranspose(matrixB);

    testSameMatrix ( testMatrixB, matrixA );

    //Cyclic Distribution
    matrixA = testMatrixA;
    matrixB = testMatrixB;
   
    IndexType chunkSize  = 1;
    DistributionPtr rowsCyclic( new CyclicDistribution( matrixA.getNumRows(), chunkSize, comm ) );
    DistributionPtr columnsCyclic( new CyclicDistribution( matrixA.getNumColumns(), chunkSize, comm ) );
    matrixA.redistribute(rowsCyclic,columnsCyclic);

    DistributionPtr rowCyclic( new CyclicDistribution( matrixB.getNumRows(), chunkSize, comm ) );
    DistributionPtr columnCyclic( new CyclicDistribution( matrixB.getNumColumns(), chunkSize, comm ) );
    matrixB.redistribute(rowCyclic,columnCyclic);

    //Transpose
    matrixA.assignTranspose(matrixA);
    matrixB.assignTranspose(matrixB);
    matrixA.assignTranspose(matrixB);

    testSameMatrix ( testMatrixB, matrixA );

    //General Block Distribution
/*    
    IndexType numProcesses = 4;

    std::vector<IndexType> localRowsSizes;
    IndexType rowsSizes[] = { 1, 1, 3, 1};
    localRowsSizes.assign( rowsSizes, rowsSizes + numProcesses );
    DistributionPtr rowsGenBlock( new GenBlockDistribution( matrixA.getNumRows(), localRowsSizes, comm ) );
    std::vector<IndexType> localColumnsSizes;
    IndexType columnsSizes[] = { 1, 1, 1, 1};
    localColumnsSizes.assign( columnsSizes, columnsSizes + numProcesses );
    DistributionPtr columnsGenBlock( new GenBlockDistribution( matrixA.getNumColumns(), localColumnsSizes, comm ) );
    matrixA.redistribute(rowsGenBlock,columnsGenBlock);

    std::vector<IndexType> localRowSizes;
    IndexType rowSizes[] = { 2, 2, 1, 1};
    localRowSizes.assign( rowSizes, rowSizes + numProcesses );
    DistributionPtr rowGenBlock( new GenBlockDistribution( matrixB.getNumRows(), localRowSizes, comm ) );
    std::vector<IndexType> localColumnSizes;
    IndexType columnSizes[] = { 1, 3, 1, 1};
    localColumnSizes.assign( columnSizes, columnSizes + numProcesses );
    DistributionPtr columnGenBlock( new GenBlockDistribution( matrixB.getNumColumns(), localColumnSizes, comm ) );
    matrixB.redistribute(rowGenBlock,columnGenBlock);
*/

    //GeneralDistribution
    // TODO: need a more general test method for working with an arbitrary number of processes
    // 4 processes
/*    IndexType distRows[] = { 0,1,3,0,2,2};
    std::vector<IndexType> row2partRows;
    row2partRows.assign( distRows, distRows + matrixA.getNumRows() );
    DistributionPtr rowsGen( new GeneralDistribution( row2partRows, matrixA.getNumRows(), comm ) );
    IndexType distColumns[] = {0,1,0,2};
    std::vector<IndexType> row2partColumns;
    row2partColumns.assign( distColumns, distColumns + matrixA.getNumColumns() );
    DistributionPtr columnsGen( new GeneralDistribution( row2partColumns, matrixA.getNumColumns(), comm ) );
    matrixA.redistribute(rowsGen,columnsGen);

    IndexType distRow[] = { 0,1,3,0,2,2};
    std::vector<IndexType> row2partRow;
    row2partRow.assign( distRow, distRow + matrixB.getNumRows() );
    DistributionPtr rowGen( new GeneralDistribution( row2partRow, matrixB.getNumRows(), comm ) );
    IndexType distColumn[] = {0,1,0,2,0,3};
    std::vector<IndexType> row2partColumn;
    row2partColumn.assign( distColumn, distColumn + matrixB.getNumColumns() );
    DistributionPtr columnGen( new GeneralDistribution( row2partColumn, matrixB.getNumColumns(), comm ) );
    matrixB.redistribute(rowGen,columnGen);


    //Transpose
    matrixA.assignTranspose(matrixA);
    matrixB.assignTranspose(matrixB);
    matrixA.assignTranspose(matrixB);
    */


    //Printings
/*
    std::cout<<"matrix A\n";
    for ( int i = 0; i < matrixA.getNumRows(); i++ ){
        for ( int j = 0; j < matrixA.getNumColumns(); j++ ){
           std::cout<< matrixA.getValue(i,j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n\n";

    std::cout<<"matrix B\n";
    for ( int i = 0; i < matrixB.getNumRows(); i++ ){
        for ( int j = 0; j < matrixB.getNumColumns(); j++ ){
           std::cout<< matrixB.getValue(i,j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n\n";
*/

}

} /* end namespace DenseMatrixTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseMatrixTest" )

LAMA_AUTO_TEST_CASE_TL( assignmentMultiplicationTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_TL( xGEMMOperationTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( cTorTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( getLocalValuesTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( getMemoryUsageTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( getNumValuesTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( getRowTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_CT( gemmTest, DenseMatrixTest, scai::lama )
LAMA_AUTO_TEST_CASE_T( matrixAddMatrixTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( maxNormTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( scaleTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( setDiagonalTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( setIdentityTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( swapTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( transposeTest, DenseMatrixTest )
LAMA_AUTO_TEST_CASE_T( typeNameTest, DenseMatrixTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
