/**
 * @file SparseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class SparseMatrixTest
 * @author Alexander Büchel
 * @date 14.03.2012
 * @since 1.0.0
 */

#include "SparseMatrixTest.hpp"

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/Scalar.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <test/TestSparseMatrices.hpp>
#include <test/Configuration.hpp>

#include <lama/NoCommunicator.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>

#include <test/TestMacros.hpp>
#include <test/SparseMatrixHelper.hpp>
#include <test/SameMatrixHelper.hpp>

using namespace lama;

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename MatrixType>, SparseMatrixTest<MatrixType>::logger,
                              "Test.SparseMatrixTest" );

/* ----------------------------------------------------------------------------- */

template<typename MatrixType>
void SparseMatrixTest<MatrixType>::setUp()
{
    std::string prefix = Configuration::getInstance().getPath();
    LAMA_LOG_INFO( logger, "prefix = " << prefix );
    m_FormattedInputFile = prefix + "/2D_poisson_256_formatted.frm";
    m_XDRInputFile = prefix + "/2D_poisson_256_xdr.frm";
    m_TestOutputFileFor = prefix + "/test_matrix_formatted.tmp.frm";
    m_TestOutputFileUnf = prefix + "/test_matrix_unformatted.tmp.frm";
    m_TestOutputFileXDR = prefix + "/test_matrix_xdr.tmp.frm";
    m_TestOutputFileMM = prefix + "/test_matrix_mm.tmp.mtx";
}

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, TypeMatrix, clearTest )LAMA_LOG_INFO( logger, "clearTest" );

typedef TypeMatrix MatrixType;
typedef typename MatrixType::ValueType ValueType;

const IndexType n = 4;

ValueType values[] =
{   4.0, 0.0, 0.0, 0.0,
    0.0, 3.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 0.0,
    0.0, 0.0, 0.0, 1.0
};

MatrixType matrix;
matrix.setRawDenseData( n, n, values );

matrix.clear();

BOOST_CHECK_EQUAL( matrix.getNumRows(), 0 );
BOOST_CHECK_EQUAL( matrix.getNumColumns(), 0 );
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, cTorTest )

LAMA_LOG_INFO( logger, "cTorTest" );

//TODO: to P_ test?
typedef typename MatrixType::ValueType ValueType;

const IndexType n = 4;

CommunicatorPtr comm = CommunicatorFactory::get( "MPI" );

DistributionPtr bdist( new BlockDistribution( n, comm ) );
DistributionPtr cdist( new CyclicDistribution( n, 1, comm ) );
DistributionPtr rdist( new NoDistribution( n ) );

//LAMA_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ")" );
//
//MatrixType matrix0( bdist );
//
//BOOST_CHECK_EQUAL( matrix0.getNumRows() , n );
//BOOST_CHECK_EQUAL( matrix0.getNumColumns() , n );

LAMA_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", bdist = " << *bdist << ")" );

MatrixType matrix1( bdist, bdist );

BOOST_CHECK_EQUAL( matrix1.getNumRows() , n );
BOOST_CHECK_EQUAL( matrix1.getNumColumns() , n );

LAMA_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", rdist = " << *bdist << ")" );

MatrixType matrix2( bdist, rdist );

BOOST_CHECK_EQUAL( matrix2.getNumRows() , n );
BOOST_CHECK_EQUAL( matrix2.getNumColumns() , n );

MatrixType matrix3;

BOOST_CHECK_EQUAL( matrix3.getNumRows() , 0 );
BOOST_CHECK_EQUAL( matrix3.getNumColumns() , 0 );

LAMA_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", cdist = " << *cdist << ")" );

MatrixType matrix4;

matrix4.setIdentity( bdist );

LAMA_LOG_INFO( logger, "Matrix( " << matrix4 << " )" );

MatrixType matrix5 ( matrix4 );

/* ToDo: provide transpose:

 MatrixType matrix6 ( matrix5 );

 BOOST_CHECK_EQUAL( matrix5.getDistribution(), matrix6.getColDistribution() );
 */
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

template<typename mt>
template<typename mt2>
void SparseMatrixTest<mt>::testConversionsImpl()
{
    LAMA_LOG_INFO( logger, "testConversionsImpl<" << typeid(mt).name() << "," << typeid(mt2).name() << ">" );

    typedef mt2 MatrixType;
    typedef typename mt2::ValueType ValueType;

    LAMA_LOG_DEBUG( logger, "conversion: n4m4" );

    CSRSparseMatrix<ValueType> csr_n4m4TestMatrix1( TestSparseMatrices::n4m4TestMatrix1<ValueType>() );

    MatrixType n4m4TestMatrix1( csr_n4m4TestMatrix1 );

    testSameMatrix( n4m4TestMatrix1, csr_n4m4TestMatrix1 );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4Inverse" );

    CSRSparseMatrix<ValueType> csr_n4m4InverseTestMatrix1( TestSparseMatrices::n4m4InverseTestMatrix1<ValueType>() );

    MatrixType n4m4InverseTestMatrix1( csr_n4m4InverseTestMatrix1 );

    testSameMatrix( n4m4InverseTestMatrix1, csr_n4m4InverseTestMatrix1 );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4Diagonal" );

    CSRSparseMatrix<ValueType> csr_n4m4DiagonalMatrix( TestSparseMatrices::n4m4DiagonalMatrix<ValueType>() );

    MatrixType n4m4DiagonalMatrix( csr_n4m4DiagonalMatrix );

    testSameMatrix( n4m4DiagonalMatrix, csr_n4m4DiagonalMatrix );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4Symmetric" );

    CSRSparseMatrix<ValueType> csr_n4m4SymmetricMatrix( TestSparseMatrices::n4m4SymmetricMatrix<ValueType>() );

    MatrixType n4m4SymmetricMatrix( csr_n4m4SymmetricMatrix );

    testSameMatrix( n4m4SymmetricMatrix, csr_n4m4SymmetricMatrix );

    LAMA_LOG_DEBUG( logger, "conversion: n4m6NonSquare" );

    CSRSparseMatrix<ValueType> csr_n4m6NoneSquareMatrix( TestSparseMatrices::n4m6NoneSquareMatrix<ValueType>() );

    MatrixType n4m6NoneSquareMatrix( csr_n4m6NoneSquareMatrix );

    testSameMatrix( n4m6NoneSquareMatrix, csr_n4m6NoneSquareMatrix );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4A1" );

    CSRSparseMatrix<ValueType> csr_n4m4MatrixA1( TestSparseMatrices::n4m4MatrixA1<ValueType>() );

    MatrixType n4m4MatrixA1( csr_n4m4MatrixA1 );

    testSameMatrix( n4m4MatrixA1, csr_n4m4MatrixA1 );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4A2" );

    CSRSparseMatrix<ValueType> csr_n4m4MatrixA2( TestSparseMatrices::n4m4MatrixA2<ValueType>() );

    MatrixType n4m4MatrixA2( csr_n4m4MatrixA2 );

    testSameMatrix( n4m4MatrixA2, csr_n4m4MatrixA2 );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4ARes" );

    CSRSparseMatrix<ValueType> csr_n4m4MatrixARes( TestSparseMatrices::n4m4MatrixARes<ValueType>() );

    MatrixType n4m4MatrixARes( csr_n4m4MatrixARes );

    testSameMatrix( n4m4MatrixARes, csr_n4m4MatrixARes );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4B1" );

    CSRSparseMatrix<ValueType> csr_n4m4MatrixB1( TestSparseMatrices::n4m4MatrixB1<ValueType>() );

    MatrixType n4m4MatrixB1( csr_n4m4MatrixB1 );

    testSameMatrix( n4m4MatrixB1, csr_n4m4MatrixB1 );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4B2" );

    CSRSparseMatrix<ValueType> csr_n4m4MatrixB2( TestSparseMatrices::n4m4MatrixB2<ValueType>() );

    MatrixType n4m4MatrixB2( csr_n4m4MatrixB2 );

    testSameMatrix( n4m4MatrixB2, csr_n4m4MatrixB2 );

    LAMA_LOG_DEBUG( logger, "conversion: n8m4Interpol" );

    CSRSparseMatrix<ValueType> csr_n8m4Interpol( TestSparseMatrices::n8m4Interpol<ValueType>() );

    MatrixType n8m4Interpol( csr_n8m4Interpol );

    testSameMatrix( n8m4Interpol, csr_n8m4Interpol );

    LAMA_LOG_DEBUG( logger, "conversion: n8m4InterpolTranspose" );

    CSRSparseMatrix<ValueType> csr_n4m8InterpolTranspose( TestSparseMatrices::n4m8InterpolTranspose<ValueType>() );

    MatrixType n4m8InterpolTranspose( csr_n4m8InterpolTranspose );

    testSameMatrix( n4m8InterpolTranspose, csr_n4m8InterpolTranspose );

    LAMA_LOG_DEBUG( logger, "conversion: n8m8Laplace1D" );

    CSRSparseMatrix<ValueType> csr_n8m8Laplace1D( TestSparseMatrices::n8m8Laplace1D<ValueType>() );

    MatrixType n8m8Laplace1D( csr_n8m8Laplace1D );

    testSameMatrix( n8m8Laplace1D, csr_n8m8Laplace1D );

    LAMA_LOG_DEBUG( logger, "conversion: n4m4Galerkin" );

    CSRSparseMatrix<ValueType> csr_n4m4Galerkin( TestSparseMatrices::n4m4Galerkin<ValueType>() );

    MatrixType n4m4Galerkin( csr_n4m4Galerkin );

    testSameMatrix( n4m4Galerkin, csr_n4m4Galerkin );
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, testConversions );

LAMA_LOG_INFO( logger, "testConversions" );

testConversionsImpl<CSRSparseMatrix<float> >();
testConversionsImpl<CSRSparseMatrix<double> >();
testConversionsImpl<ELLSparseMatrix<float> >();
testConversionsImpl<ELLSparseMatrix<double> >();
testConversionsImpl<JDSSparseMatrix<float> >();
testConversionsImpl<JDSSparseMatrix<double> >();
testConversionsImpl<DIASparseMatrix<float> >();
testConversionsImpl<DIASparseMatrix<double> >();
testConversionsImpl<COOSparseMatrix<float> >();
testConversionsImpl<COOSparseMatrix<double> >();
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

template<typename MatrixType>
void SparseMatrixTest<MatrixType>::matrixMultTestImpl( const Matrix& a, const Matrix& b, const Matrix& result )
{
    typedef typename MatrixType::ValueType ValueType;

    LAMA_LOG_INFO( logger,
                   "matrixMultTestImpl: verify a * b = result, with a = " << a << ", b = " << b << ", result = " << result );

    MatrixType csrResult;

    MatrixType A( a );
    MatrixType B( b );

    csrResult = A * B;

    BOOST_CHECK_EQUAL( result.getNumRows(), csrResult.getNumRows() );
    BOOST_CHECK_EQUAL( result.getNumColumns(), csrResult.getNumColumns() );

    MatrixType res( result );

    for ( IndexType i = 0; i < csrResult.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < csrResult.getNumColumns(); ++j )
        {
            Scalar s1 = csrResult.getValue( i, j );
            Scalar s2 = res.getValue( i, j );
            BOOST_CHECK_CLOSE( s1.getValue<ValueType>(), s2.getValue<ValueType>(), 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, testMultiplication )

LAMA_LOG_INFO( logger, "testMultiplication" );

typedef typename MatrixType::ValueType ValueType;

//TODO: crashes, because of missing Ctor in SparseMatrix with DenseMatrix as argument
//    CSRSparseMatrix<double> randomMatrix =
//          TestSparseMatrices::n4m4TestMatrix1<double>();
//
//    IndexType n = randomMatrix.getNumRows();
//
//    DenseMatrix<double> rDM( randomMatrix );
//    DenseMatrix<double> randomDenseMatrix( randomMatrix );
//    DenseMatrixOps Dmops;
//    Dmops.invert( randomDenseMatrix );

//TODO: here

//CSRSparseMatrix<double> InvertedMatrix( randomDenseMatrix );

//    CSRSparseMatrix<double> identity =
//        TestSparseMatrices::nnIdentityMatrix<double>(n);
//    matrixEqualityCheck( result, IdentityDenseMatrix );

const CSRSparseMatrix<ValueType> matrixA =
    TestSparseMatrices::n4m4MatrixA1<ValueType>();

const CSRSparseMatrix<ValueType> matrixB =
    TestSparseMatrices::n4m4MatrixA2<ValueType>();

const CSRSparseMatrix<ValueType> matrixC =
    TestSparseMatrices::n4m4MatrixARes<ValueType>();

LAMA_LOG_INFO( logger, "verify: n4m4MatrixA1 * n4m4MatrixA2 = n4m4MatrixARes" );

matrixMultTestImpl( matrixA, matrixB, matrixC );

const CSRSparseMatrix<ValueType> matrixD =
    TestSparseMatrices::n4m4MatrixB1<ValueType>();

const CSRSparseMatrix<ValueType> matrixE =
    TestSparseMatrices::n4m4MatrixB2<ValueType>();

const CSRSparseMatrix<ValueType> matrixF =
    TestSparseMatrices::n4m4MatrixBRes<ValueType>();

LAMA_LOG_INFO( logger, "verify: n4m4MatrixB1 * n4m4MatrixB2 = n4m4MatrixBRes" );

matrixMultTestImpl( matrixD, matrixE, matrixF );

const CSRSparseMatrix<ValueType> matrixG =
    TestSparseMatrices::n4m4MatrixC1<ValueType>();

const CSRSparseMatrix<ValueType> matrixGxG =
    TestSparseMatrices::n4m4MatrixCRes<ValueType>();

LAMA_LOG_INFO( logger, "verify: n4m4MatrixC1 * n4m4MatrixC1 = n4m4MatrixCRes" );

matrixMultTestImpl( matrixG, matrixG, matrixGxG );

const CSRSparseMatrix<ValueType> matrixA1 =
    TestSparseMatrices::n6m4MatrixD1<ValueType>();

const CSRSparseMatrix<ValueType> matrixB1 =
    TestSparseMatrices::n4m6MatrixD2<ValueType>();

const CSRSparseMatrix<ValueType> matrixC1 =
    TestSparseMatrices::n6m6MatrixDRes<ValueType>();

LAMA_LOG_INFO( logger, "verify: n6m4MatrixD1 * n4m4MatrixD2 = n4m4MatrixDRes" );

matrixMultTestImpl( matrixA1, matrixB1, matrixC1 );

const CSRSparseMatrix<ValueType> matrixD1 =
    TestSparseMatrices::n6m4MatrixE1<ValueType>();

const CSRSparseMatrix<ValueType> matrixE1 =
    TestSparseMatrices::n4m3MatrixE2<ValueType>();

const CSRSparseMatrix<ValueType> matrixF1 =
    TestSparseMatrices::n6m3MatrixERes<ValueType>();

LAMA_LOG_INFO( logger, "verify: n6m4MatrixE1 * n4m3MatrixE2 = n6m3MatrixDRes" );

matrixMultTestImpl( matrixD1, matrixE1, matrixF1 );

const CSRSparseMatrix<ValueType> laplace1Dmatrix =
    TestSparseMatrices::n8m8Laplace1D<ValueType>();

CSRSparseMatrix<ValueType> interpolationMatrix( TestSparseMatrices::n8m4Interpol<ValueType>() );

CSRSparseMatrix<ValueType> interpolationMatrixTransp =
    TestSparseMatrices::n8m4Interpol<ValueType>();

interpolationMatrixTransp.assignTranspose( interpolationMatrix );

LAMA_LOG_INFO( logger, "verify: n8m8Laplace1D * n8m4Interpol_Trans = n8m4GalerkinTemp" );

matrixMultTestImpl( laplace1Dmatrix, interpolationMatrix, TestSparseMatrices::n8m4GalerkinTemp<ValueType>() );

CSRSparseMatrix<ValueType> laplaceTimesInterpol;
laplaceTimesInterpol = laplace1Dmatrix * interpolationMatrix;

CSRSparseMatrix<ValueType> coarseGridGalerkin;
coarseGridGalerkin = interpolationMatrixTransp * laplaceTimesInterpol;

matrixMultTestImpl( interpolationMatrixTransp, laplaceTimesInterpol, TestSparseMatrices::n4m4Galerkin<ValueType>() );
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

template<typename MatrixType>
void SparseMatrixTest<MatrixType>::matrixEqualityCheck( const MatrixType& a, const MatrixType& b )
{
    typedef typename MatrixType::ValueType ValueType;

    BOOST_CHECK_EQUAL( a.getNumRows(), b.getNumRows() );
    BOOST_CHECK_EQUAL( a.getNumColumns(), b.getNumColumns() );

    for ( IndexType i = 0; i < a.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < b.getNumColumns(); ++j )
        {
            Scalar s1 = a.getValue( i, j );
            Scalar s2 = b.getValue( i, j );
            LAMA_CHECK_SCALAR_CLOSE( s1, s2, ValueType, 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, MatrixExpressionTest )

LAMA_LOG_INFO( logger, "MatrixExpressionTest" );

//Creating test objects
typedef typename MatrixType::ValueType ValueType;

MatrixType testmatrix;
const IndexType n = 4;

MatrixType matrixA( TestSparseMatrices::n4m4MatrixB1<ValueType>() );
MatrixType matrixB( TestSparseMatrices::n4m4MatrixB2<ValueType>() );
MatrixType matrixC( TestSparseMatrices::n4m4MatrixB1<ValueType>() );

Scalar s = 2.0;
Scalar t = 3.0;

//Expression-test A*B
MatrixType Result1( TestSparseMatrices::n4m4MatrixBRes<ValueType>() );
testmatrix = matrixA * matrixB;
matrixEqualityCheck( testmatrix, Result1 );

//Expression-test a*A*B
ValueType valuesResult2[] =
{   8.0, 0.0, 0.0, 0.0,
    0.0,12.0, 0.0, 0.0,
    0.0, 0.0,12.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};

MatrixType Result2;
Result2.setRawDenseData( n, n, valuesResult2 );
testmatrix.clear();
testmatrix = s * matrixA * matrixB;
matrixEqualityCheck( testmatrix, Result2 );

//Expression-test A*a*B
testmatrix.clear();
testmatrix = matrixA * s * matrixB;
matrixEqualityCheck( testmatrix, Result2 );

//Expression-test A*B*a
testmatrix.clear();
testmatrix = matrixA * matrixB * s;
matrixEqualityCheck( testmatrix, Result2 );

//Expression b*C
ValueType valuesResult3[] =
{   2.0, 0.0, 0.0, 0.0,
    0.0, 4.0, 0.0, 0.0,
    0.0, 0.0, 6.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};

MatrixType result3;
result3.setRawDenseData( n, n, valuesResult3 );

testmatrix.clear();
testmatrix = s * matrixA;
matrixEqualityCheck( testmatrix, result3 );

//Expression C*b
testmatrix.clear();
testmatrix = matrixA * s;
matrixEqualityCheck( testmatrix, result3 );

//Expression 1*C
testmatrix.clear();
Scalar u = 1.0;
testmatrix = u * matrixA;
matrixEqualityCheck( testmatrix, matrixA );

//Expression C*1
testmatrix.clear();
testmatrix = matrixA * u;
matrixEqualityCheck( testmatrix, matrixA );
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, MatrixCtorExpressionTest )

LAMA_LOG_INFO( logger, "MatrixCtorExpressionTest" );

//Creating test objects
typedef typename MatrixType::ValueType ValueType;

const IndexType n = 4;

MatrixType matrixA( TestSparseMatrices::n4m4MatrixB1<ValueType>() );
MatrixType matrixB( TestSparseMatrices::n4m4MatrixB2<ValueType>() );
MatrixType matrixC( TestSparseMatrices::n4m4MatrixB1<ValueType>() );

Scalar s = 2.0;
Scalar t = 3.0;
Scalar u = 1.0;

//Expression-test A*B
MatrixType Result1( TestSparseMatrices::n4m4MatrixBRes<ValueType>() );
MatrixType testmatrix1 ( matrixA * matrixB );
matrixEqualityCheck( testmatrix1, Result1 );

//Expression b*C
ValueType valuesResult3[] =
{   2.0, 0.0, 0.0, 0.0,
    0.0, 4.0, 0.0, 0.0,
    0.0, 0.0, 6.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};
MatrixType result3;
result3.setRawDenseData( n, n, valuesResult3 );
MatrixType testmatrix2( s * matrixA );
matrixEqualityCheck( testmatrix2, result3 );

//Expression C*b
MatrixType testmatrix3( matrixA * s );
matrixEqualityCheck( testmatrix3, result3 );

//Expression 1*B
MatrixType testmatrix7( u * matrixA );
matrixEqualityCheck( testmatrix7, matrixA );

//Expression B*1
MatrixType testmatrix8( matrixA * u );
matrixEqualityCheck( testmatrix8, matrixA );

//Expression-test a*A*B
ValueType valuesResult2[] =
{   8.0, 0.0, 0.0, 0.0,
    0.0,12.0, 0.0, 0.0,
    0.0, 0.0,12.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};

MatrixType result2;
result2.setRawDenseData( n, n, valuesResult2 );

MatrixType testmatrix4 ( s * matrixA * matrixB );
matrixEqualityCheck( testmatrix4, result2 );

//Expression-test A*a*B
MatrixType testmatrix5 ( matrixA * s * matrixB );
matrixEqualityCheck( testmatrix5, result2 );

//Expression-test A*B*a
MatrixType testmatrix6 ( matrixA * matrixB * s );
matrixEqualityCheck( testmatrix6, result2 );

LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, writeAtTest )

LAMA_WRITEAT_TEST( mMatrix );

LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, scaleTest )MatrixType matrixA( TestSparseMatrices::n4m4MatrixA1<double>() );

double valuesResult[] =
{   1.2f, 0.0f, 0.0f, 0.8f,
    1.4f, 0.8f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.8f, 0.8f,
    0.4f, 1.0f, 0.0f, 0.6f
};

MatrixType result;
result.setRawDenseData( 4, 4, valuesResult );

Scalar s = 2.0;

matrixA.scale( s );

matrixEqualityCheck( matrixA, result );
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( SparseMatrixTest ) {
    clearTest();
    cTorTest();
    testConversions();
    testMultiplication();
    MatrixExpressionTest();
    MatrixCtorExpressionTest();
    writeAtTest();
    scaleTest();
}

/* ----------------------------------------------------------------------------- */

template class SparseMatrixTest<CSRSparseMatrix<double> > ;
template class SparseMatrixTest<CSRSparseMatrix<float> > ;
template class SparseMatrixTest<ELLSparseMatrix<double> > ;
template class SparseMatrixTest<ELLSparseMatrix<float> > ;
template class SparseMatrixTest<JDSSparseMatrix<double> > ;
template class SparseMatrixTest<JDSSparseMatrix<float> > ;
template class SparseMatrixTest<COOSparseMatrix<double> > ;
template class SparseMatrixTest<COOSparseMatrix<float> > ;
template class SparseMatrixTest<DIASparseMatrix<double> > ;
template class SparseMatrixTest<DIASparseMatrix<float> > ;
template class SparseMatrixTest<DenseMatrix<double> > ;
template class SparseMatrixTest<DenseMatrix<float> > ;
