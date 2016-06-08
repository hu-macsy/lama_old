/**
 * @file SparseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains the implementation of the class SparseMatrixTest
 * @author Alexander Büchel
 * @date 14.03.2012
 */

#include <scai/lama/test/SparseMatrixTest.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/test/TestSparseMatrices.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/SparseMatrixHelper.hpp>
#include <scai/lama/test/SameMatrixHelper.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename MatrixType>, SparseMatrixTest<MatrixType>::logger,
                              "Test.SparseMatrixTest" )

/* ----------------------------------------------------------------------------- */

template<typename MatrixType>
void SparseMatrixTest<MatrixType>::setUp()
{
    std::string prefix = scai::test::Configuration::getPath();
    SCAI_LOG_INFO( logger, "prefix = " << prefix );
    m_FormattedInputFile = prefix + "/2D_poisson_256_formatted.frm";
    m_XDRInputFile = prefix + "/2D_poisson_256_xdr.frm";
    m_TestOutputFileFor = prefix + "/test_matrix_formatted.tmp.frm";
    m_TestOutputFileUnf = prefix + "/test_matrix_unformatted.tmp.frm";
    m_TestOutputFileXDR = prefix + "/test_matrix_xdr.tmp.frm";
    m_TestOutputFileMM = prefix + "/test_matrix_mm.tmp.mtx";
}

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, clearTest )

SCAI_LOG_INFO( logger, "clearTest" )

typedef typename MatrixType::MatrixValueType ValueType;

const IndexType n = 4;

ValueType values[] =
{
    4.0, 0.0, 0.0, 0.0,
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

SCAI_LOG_INFO( logger, "cTorTest" );

//TODO: to P_ test?
const IndexType n = 4;

CommunicatorPtr comm = Communicator::getCommunicatorPtr( Communicator::MPI );

DistributionPtr bdist( new BlockDistribution( n, comm ) );
DistributionPtr cdist( new CyclicDistribution( n, 1, comm ) );
DistributionPtr rdist( new NoDistribution( n ) );

//SCAI_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ")" );
//
//MatrixType matrix0( bdist );
//
//BOOST_CHECK_EQUAL( matrix0.getNumRows() , n );
//BOOST_CHECK_EQUAL( matrix0.getNumColumns() , n );

SCAI_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", bdist = " << *bdist << ")" );

MatrixType matrix1( bdist, bdist );

BOOST_CHECK_EQUAL( matrix1.getNumRows() , n );
BOOST_CHECK_EQUAL( matrix1.getNumColumns() , n );

SCAI_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", rdist = " << *bdist << ")" );

MatrixType matrix2( bdist, rdist );

BOOST_CHECK_EQUAL( matrix2.getNumRows() , n );
BOOST_CHECK_EQUAL( matrix2.getNumColumns() , n );

MatrixType matrix3;

BOOST_CHECK_EQUAL( matrix3.getNumRows() , 0 );
BOOST_CHECK_EQUAL( matrix3.getNumColumns() , 0 );

SCAI_LOG_INFO( logger, "Matrix( bdist = " << *bdist << ", cdist = " << *cdist << ")" );

MatrixType matrix4;

matrix4.setIdentity( bdist );

SCAI_LOG_INFO( logger, "Matrix( " << matrix4 << " )" );

MatrixType matrix5 ( matrix4 );

/* ToDo: provide transpose:

 MatrixType matrix6 ( matrix5 );

 BOOST_CHECK_EQUAL( matrix5.getDistribution(), matrix6.getColDistribution() );
 */
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

template<typename MatrixType1>
template<typename MatrixType2>
void SparseMatrixTest<MatrixType1>::testConversionsImpl()
{
    SCAI_LOG_INFO( logger, "testConversionsImpl<" << typeid( MatrixType1 ).name() << "," << typeid( MatrixType2 ).name() << ">" );
    typedef MatrixType2 MatrixType;
    typedef typename MatrixType2::MatrixValueType ValueType;
    SCAI_LOG_DEBUG( logger, "conversion: n4m4" );
    CSRSparseMatrix<ValueType> csr_n4m4TestMatrix1( TestSparseMatrices::n4m4TestMatrix1<ValueType>() );
    MatrixType n4m4TestMatrix1( csr_n4m4TestMatrix1 );
    testSameMatrix( n4m4TestMatrix1, csr_n4m4TestMatrix1 );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4Inverse" );
    CSRSparseMatrix<ValueType> csr_n4m4InverseTestMatrix1( TestSparseMatrices::n4m4InverseTestMatrix1<ValueType>() );
    MatrixType n4m4InverseTestMatrix1( csr_n4m4InverseTestMatrix1 );
    testSameMatrix( n4m4InverseTestMatrix1, csr_n4m4InverseTestMatrix1 );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4Diagonal" );
    CSRSparseMatrix<ValueType> csr_n4m4DiagonalMatrix( TestSparseMatrices::n4m4DiagonalMatrix<ValueType>() );
    MatrixType n4m4DiagonalMatrix( csr_n4m4DiagonalMatrix );
    testSameMatrix( n4m4DiagonalMatrix, csr_n4m4DiagonalMatrix );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4Symmetric" );
    CSRSparseMatrix<ValueType> csr_n4m4SymmetricMatrix( TestSparseMatrices::n4m4SymmetricMatrix<ValueType>() );
    MatrixType n4m4SymmetricMatrix( csr_n4m4SymmetricMatrix );
    testSameMatrix( n4m4SymmetricMatrix, csr_n4m4SymmetricMatrix );
    SCAI_LOG_DEBUG( logger, "conversion: n4m6NonSquare" );
    CSRSparseMatrix<ValueType> csr_n4m6NoneSquareMatrix( TestSparseMatrices::n4m6NoneSquareMatrix<ValueType>() );
    MatrixType n4m6NoneSquareMatrix( csr_n4m6NoneSquareMatrix );
    testSameMatrix( n4m6NoneSquareMatrix, csr_n4m6NoneSquareMatrix );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4A1" );
    CSRSparseMatrix<ValueType> csr_n4m4MatrixA1( TestSparseMatrices::n4m4MatrixA1<ValueType>() );
    MatrixType n4m4MatrixA1( csr_n4m4MatrixA1 );
    testSameMatrix( n4m4MatrixA1, csr_n4m4MatrixA1 );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4A2" );
    CSRSparseMatrix<ValueType> csr_n4m4MatrixA2( TestSparseMatrices::n4m4MatrixA2<ValueType>() );
    MatrixType n4m4MatrixA2( csr_n4m4MatrixA2 );
    testSameMatrix( n4m4MatrixA2, csr_n4m4MatrixA2 );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4ARes" );
    CSRSparseMatrix<ValueType> csr_n4m4MatrixARes( TestSparseMatrices::n4m4MatrixARes<ValueType>() );
    MatrixType n4m4MatrixARes( csr_n4m4MatrixARes );
    testSameMatrix( n4m4MatrixARes, csr_n4m4MatrixARes );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4B1" );
    CSRSparseMatrix<ValueType> csr_n4m4MatrixB1( TestSparseMatrices::n4m4MatrixB1<ValueType>() );
    MatrixType n4m4MatrixB1( csr_n4m4MatrixB1 );
    testSameMatrix( n4m4MatrixB1, csr_n4m4MatrixB1 );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4B2" );
    CSRSparseMatrix<ValueType> csr_n4m4MatrixB2( TestSparseMatrices::n4m4MatrixB2<ValueType>() );
    MatrixType n4m4MatrixB2( csr_n4m4MatrixB2 );
    testSameMatrix( n4m4MatrixB2, csr_n4m4MatrixB2 );
    SCAI_LOG_DEBUG( logger, "conversion: n8m4Interpol" );
    CSRSparseMatrix<ValueType> csr_n8m4Interpol( TestSparseMatrices::n8m4Interpol<ValueType>() );
    MatrixType n8m4Interpol( csr_n8m4Interpol );
    testSameMatrix( n8m4Interpol, csr_n8m4Interpol );
    SCAI_LOG_DEBUG( logger, "conversion: n8m4InterpolTranspose" );
    CSRSparseMatrix<ValueType> csr_n4m8InterpolTranspose( TestSparseMatrices::n4m8InterpolTranspose<ValueType>() );
    MatrixType n4m8InterpolTranspose( csr_n4m8InterpolTranspose );
    testSameMatrix( n4m8InterpolTranspose, csr_n4m8InterpolTranspose );
    SCAI_LOG_DEBUG( logger, "conversion: n8m8Laplace1D" );
    CSRSparseMatrix<ValueType> csr_n8m8Laplace1D( TestSparseMatrices::n8m8Laplace1D<ValueType>() );
    MatrixType n8m8Laplace1D( csr_n8m8Laplace1D );
    testSameMatrix( n8m8Laplace1D, csr_n8m8Laplace1D );
    SCAI_LOG_DEBUG( logger, "conversion: n4m4Galerkin" );
    CSRSparseMatrix<ValueType> csr_n4m4Galerkin( TestSparseMatrices::n4m4Galerkin<ValueType>() );
    MatrixType n4m4Galerkin( csr_n4m4Galerkin );
    testSameMatrix( n4m4Galerkin, csr_n4m4Galerkin );
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, testConversions )

SCAI_LOG_INFO( logger, "testConversions" );

typedef SCAI_TEST_TYPE ValueType;

testConversionsImpl<CSRSparseMatrix<ValueType> >();
testConversionsImpl<ELLSparseMatrix<ValueType> >();
testConversionsImpl<JDSSparseMatrix<ValueType> >();
testConversionsImpl<DIASparseMatrix<ValueType> >();
testConversionsImpl<COOSparseMatrix<ValueType> >();

LAMA_COMMON_TEST_CASE_TEMPLATE_END();

template<typename MatrixType>
void SparseMatrixTest<MatrixType>::matrixMultTestImpl( const Matrix& a, const Matrix& b, const Matrix& result )
{
    SCAI_LOG_INFO( logger,
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
            SCAI_CHECK_CLOSE( csrResult.getValue( i, j ), res.getValue( i, j ), 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, testMultiplication )

SCAI_LOG_INFO( logger, "testMultiplication" );

typedef typename MatrixType::MatrixValueType ValueType;

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

SCAI_LOG_INFO( logger, "verify: n4m4MatrixA1 * n4m4MatrixA2 = n4m4MatrixARes" );

matrixMultTestImpl( matrixA, matrixB, matrixC );

const CSRSparseMatrix<ValueType> matrixD =
    TestSparseMatrices::n4m4MatrixB1<ValueType>();

const CSRSparseMatrix<ValueType> matrixE =
    TestSparseMatrices::n4m4MatrixB2<ValueType>();

const CSRSparseMatrix<ValueType> matrixF =
    TestSparseMatrices::n4m4MatrixBRes<ValueType>();

SCAI_LOG_INFO( logger, "verify: n4m4MatrixB1 * n4m4MatrixB2 = n4m4MatrixBRes" );

matrixMultTestImpl( matrixD, matrixE, matrixF );

const CSRSparseMatrix<ValueType> matrixG =
    TestSparseMatrices::n4m4MatrixC1<ValueType>();

const CSRSparseMatrix<ValueType> matrixGxG =
    TestSparseMatrices::n4m4MatrixCRes<ValueType>();

SCAI_LOG_INFO( logger, "verify: n4m4MatrixC1 * n4m4MatrixC1 = n4m4MatrixCRes" );

matrixMultTestImpl( matrixG, matrixG, matrixGxG );

const CSRSparseMatrix<ValueType> matrixA1 =
    TestSparseMatrices::n6m4MatrixD1<ValueType>();

const CSRSparseMatrix<ValueType> matrixB1 =
    TestSparseMatrices::n4m6MatrixD2<ValueType>();

const CSRSparseMatrix<ValueType> matrixC1 =
    TestSparseMatrices::n6m6MatrixDRes<ValueType>();

SCAI_LOG_INFO( logger, "verify: n6m4MatrixD1 * n4m4MatrixD2 = n4m4MatrixDRes" );

matrixMultTestImpl( matrixA1, matrixB1, matrixC1 );

const CSRSparseMatrix<ValueType> matrixD1 =
    TestSparseMatrices::n6m4MatrixE1<ValueType>();

const CSRSparseMatrix<ValueType> matrixE1 =
    TestSparseMatrices::n4m3MatrixE2<ValueType>();

const CSRSparseMatrix<ValueType> matrixF1 =
    TestSparseMatrices::n6m3MatrixERes<ValueType>();

SCAI_LOG_INFO( logger, "verify: n6m4MatrixE1 * n4m3MatrixE2 = n6m3MatrixDRes" );

matrixMultTestImpl( matrixD1, matrixE1, matrixF1 );

const CSRSparseMatrix<ValueType> laplace1Dmatrix =
    TestSparseMatrices::n8m8Laplace1D<ValueType>();

CSRSparseMatrix<ValueType> interpolationMatrix( TestSparseMatrices::n8m4Interpol<ValueType>() );

CSRSparseMatrix<ValueType> interpolationMatrixTransp =
    TestSparseMatrices::n8m4Interpol<ValueType>();

interpolationMatrixTransp.assignTranspose( interpolationMatrix );

SCAI_LOG_INFO( logger, "verify: n8m8Laplace1D * n8m4Interpol_Trans = n8m4GalerkinTemp" );

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
    BOOST_CHECK_EQUAL( a.getNumRows(), b.getNumRows() );
    BOOST_CHECK_EQUAL( a.getNumColumns(), b.getNumColumns() );

    for ( IndexType i = 0; i < a.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < b.getNumColumns(); ++j )
        {
            SCAI_CHECK_CLOSE( a.getValue( i, j ), b.getValue( i, j ), 1 );
        }
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, MatrixExpressionTest )

SCAI_LOG_INFO( logger, "MatrixExpressionTest" )

//Creating test objects
typedef typename MatrixType::MatrixValueType ValueType;

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
{
    8.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 0.0, 0.0,
    0.0, 0.0, 12.0, 0.0,
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
{
    2.0, 0.0, 0.0, 0.0,
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

SCAI_LOG_INFO( logger, "MatrixCtorExpressionTest" );

//Creating test objects
typedef typename MatrixType::MatrixValueType ValueType;

const IndexType n = 4;

MatrixType matrixA( TestSparseMatrices::n4m4MatrixB1<ValueType>() );
MatrixType matrixB( TestSparseMatrices::n4m4MatrixB2<ValueType>() );
MatrixType matrixC( TestSparseMatrices::n4m4MatrixB1<ValueType>() );

Scalar s = 2.0;
Scalar t = 3.0;
Scalar u = 1.0;

//Expression-test A*B
MatrixType Result1( TestSparseMatrices::n4m4MatrixBRes<ValueType>() );
MatrixType testmatrix1 ( matrixA* matrixB );
matrixEqualityCheck( testmatrix1, Result1 );

//Expression b*C
ValueType valuesResult3[] =
{
    2.0, 0.0, 0.0, 0.0,
    0.0, 4.0, 0.0, 0.0,
    0.0, 0.0, 6.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};
MatrixType result3;
result3.setRawDenseData( n, n, valuesResult3 );
MatrixType testmatrix2( s* matrixA );
matrixEqualityCheck( testmatrix2, result3 );

//Expression C*b
MatrixType testmatrix3( matrixA* s );
matrixEqualityCheck( testmatrix3, result3 );

//Expression 1*B
MatrixType testmatrix7( u* matrixA );
matrixEqualityCheck( testmatrix7, matrixA );

//Expression B*1
MatrixType testmatrix8( matrixA* u );
matrixEqualityCheck( testmatrix8, matrixA );

//Expression-test a*A*B
ValueType valuesResult2[] =
{
    8.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 0.0, 0.0,
    0.0, 0.0, 12.0, 0.0,
    0.0, 0.0, 0.0, 8.0
};

MatrixType result2;
result2.setRawDenseData( n, n, valuesResult2 );

MatrixType testmatrix4 ( s* matrixA* matrixB );
matrixEqualityCheck( testmatrix4, result2 );

//Expression-test A*a*B
MatrixType testmatrix5 ( matrixA* s* matrixB );
matrixEqualityCheck( testmatrix5, result2 );

//Expression-test A*B*a
MatrixType testmatrix6 ( matrixA* matrixB* s );
matrixEqualityCheck( testmatrix6, result2 );

LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, writeAtTest )
SCAI_COMMON_WRITEAT_TEST( mMatrix );
LAMA_COMMON_TEST_CASE_TEMPLATE_END()

/* ----------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( SparseMatrixTest, MatrixType, scaleTest )

typedef SCAI_TEST_TYPE ValueType;

MatrixType matrixA( TestSparseMatrices::n4m4MatrixA1<ValueType>() );
ValueType valuesResult[] =
{
    1.2f, 0.0f, 0.0f, 0.8f,
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

template<typename StorageType>                                                                                     \
void SparseMatrixTest<StorageType>::runTests()
{
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

SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, COOSparseMatrix, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, CSRSparseMatrix, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, DIASparseMatrix, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, ELLSparseMatrix, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, JDSSparseMatrix, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS_II( SparseMatrixTest, DenseMatrix, SCAI_ARITHMETIC_HOST )

