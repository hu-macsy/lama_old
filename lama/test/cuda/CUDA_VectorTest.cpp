/**
 * @file CUDA_VectorTest.cpp
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
 * @brief Contains the implementation of the class CUDA_VectorTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 01.03.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/HostReadAccess.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/DenseVector.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/matrix/ELLSparseMatrix.hpp>
#include <lama/matrix/COOSparseMatrix.hpp>
#include <lama/matrix/DIASparseMatrix.hpp>
#include <lama/matrix/JDSSparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

#include <lama/cuda/CUDAHostContextManager.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/Distribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>

#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>

#include <test/Configuration.hpp>
#include <test/TestSparseMatrices.hpp>
#include <test/EquationHelper.hpp>
#include <test/cuda/CUDAContext.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <typeinfo>

using namespace lama;
using namespace memory;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct CUDA_VectorTestConfig
{
    CUDA_VectorTestConfig()
    {
        comm = CommunicatorFactory::get();
        std::srand( time( NULL ) );
    }

    ~CUDA_VectorTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( CUDA_VectorTest, CUDA_VectorTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.CUDA_VectorTest" )

/* ------------------------------------------------------------------------- */

void doMatrixTimesVector( Vector& y, const Matrix& A, const Vector& x, const Vector& corResult )
{
    y = 0.0;
    y = A * x;

    for ( IndexType i = 0; i < corResult.size(); ++i )
    {
        BOOST_CHECK_EQUAL( corResult.getValue( i ), y.getValue( i ) );
    }

    y = 0.0;
}

void doMatrixTimesVectorSyncAsyncTests( Vector& y, Matrix& A, const Vector& x, const Vector& corResult )
{
    Matrix::SyncKind saveSyncKind = A.getCommunicationKind();
    //1. Synchronous
    A.setCommunicationKind( Matrix::SYNCHRONOUS );
    LAMA_LOG_INFO( logger, "Communicate sync" )
    doMatrixTimesVector( y, A, x, corResult );
    //2. Asynchronous
    A.setCommunicationKind( Matrix::ASYNCHRONOUS );
    LAMA_LOG_INFO( logger, "Communicate async" )
    doMatrixTimesVector( y, A, x, corResult );
    //reset to original value
    A.setCommunicationKind( saveSyncKind );
}

template<typename MatrixType>
void doMatrixTimesVectorLocationTests( Vector& y, MatrixType& A, const Vector& x, const Vector& corResult )
{
    //1. Host, Host
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    LAMA_LOG_INFO( logger, "Run local on Host, halo on Host" )
    A.setContext( hostContext, hostContext );
    doMatrixTimesVectorSyncAsyncTests( y, A, x, corResult );
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    //2. CUDA, Host
    LAMA_LOG_INFO( logger, "Run local on CUDA, halo on Host" )
    A.setContext( cudaContext, hostContext );
    doMatrixTimesVectorSyncAsyncTests( y, A, x, corResult );
    //3. Host, CUDA
    LAMA_LOG_INFO( logger, "Run local on Host, halo on Cuda" )
    A.setContext( hostContext, cudaContext );
    doMatrixTimesVectorSyncAsyncTests( y, A, x, corResult );
    //4. CUDA, CUDA
    LAMA_LOG_INFO( logger, "Run local on CUDA, halo on CUDA" )
    A.setContext( cudaContext, cudaContext );
    doMatrixTimesVectorSyncAsyncTests( y, A, x, corResult );
    //reset to defaults
    A.setContext( hostContext, hostContext );
}

template<typename MatrixType>
void matrixTimesVectorTestImpl()
{
    typedef typename MatrixType::MatrixValueType ValueType;
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    CUDAHostContextManager::setAsCurrent( cuda );
    LAMA_LOG_INFO( logger, "set CUDAHostContext as CUDAHostContextManager" )
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    shared_ptr<Distribution> dist( new BlockDistribution( vectorSize, comm ) );
    MatrixType matrixTypeMatrix( dist, dist );
    DenseVector<ValueType> denseVector( dist, 1.0 );
    DenseVector<ValueType> denseTemp( dist, 0.0 );
    // Run matrix-vector multiplication with zero matrix
    doMatrixTimesVectorLocationTests( denseTemp, matrixTypeMatrix, denseVector, denseTemp );
    int numRows = 4 * size;
    int numCols = 4 * size;
    DenseVector<ValueType> denseCorrectResult2( dist, 0.0 );
    LAMAArray<ValueType>& localDenseCorrectResult2 = denseCorrectResult2.getLocalValues();
    scoped_array<ValueType> values( new ValueType[numRows * numCols] );
    {
        HostWriteAccess<ValueType> localDenseCorrectResult2Access( localDenseCorrectResult2 );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType rowSum = 0.0;

            for ( IndexType j = 0; j < numCols; ++j )
            {
                ValueType value = 0.0;

                if ( j == i || j + size == i || j - size == i || j + 2 * size == i || j - 2 * size == i
                        || j + ( numRows - 1 ) == i || j - ( numRows - 1 ) == i )
                {
                    value = 1000.0f * ( i + 1 ) + ( j + 1 );
                }

                values[i * numCols + j] = value;
                rowSum += value;
            }

            if ( dist->isLocal( i ) )
            {
                localDenseCorrectResult2Access[dist->global2local( i )] = rowSum;
            }
        }
    }
    MatrixType repM;
    repM.setRawDenseData( numRows, numCols, values.get() );
    DenseVector<ValueType> denseVector0( vectorSize, 1.0 );
    DenseVector<ValueType> denseResult0( vectorSize, 0.0 );
    doMatrixTimesVectorLocationTests( denseResult0, repM, denseVector0, denseCorrectResult2 );
    MatrixType matrixTypeMatrix2( repM, dist, dist );
    doMatrixTimesVectorLocationTests( denseTemp, matrixTypeMatrix2, denseVector, denseCorrectResult2 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixTimesVectorTest, ValueType, test_types )
{
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for CSR" )
    matrixTimesVectorTestImpl< CSRSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for ELL" )
    matrixTimesVectorTestImpl< ELLSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for JDS" )
    matrixTimesVectorTestImpl< JDSSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for COO" )
    matrixTimesVectorTestImpl< COOSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for DIA" )
    matrixTimesVectorTestImpl< DIASparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "matrixTimesVectorTest: call implementation for Dense" )
    matrixTimesVectorTestImpl< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */

void doVectorTimesMatrix( Vector& y, const Matrix& A, const Vector& x, const Vector& corResult )
{
    y = 0.0;
    y = x * A;

    for ( IndexType i = 0; i < corResult.size(); ++i )
    {
        BOOST_CHECK_EQUAL( corResult.getValue( i ), y.getValue( i ) );
    }

    y = 0.0;
}

void doVectorTimesMatrixSyncAsyncTests( Vector& y, Matrix& A, const Vector& x, const Vector& corResult )
{
    Matrix::SyncKind saveSyncKind = A.getCommunicationKind();
    //1. Synchronous
    A.setCommunicationKind( Matrix::SYNCHRONOUS );
    LAMA_LOG_INFO( logger, "Communicate sync" )
    doVectorTimesMatrix( y, A, x, corResult );
    //2. Asynchronous
    A.setCommunicationKind( Matrix::ASYNCHRONOUS );
    LAMA_LOG_INFO( logger, "Communicate async" )
    doVectorTimesMatrix( y, A, x, corResult );
    //reset to original value
    A.setCommunicationKind( saveSyncKind );
}

template<typename MatrixType>
void doVectorTimesMatrixLocationTests( Vector& y, MatrixType& A, const Vector& x, const Vector& corResult )
{
    //1. Host, Host
    ContextPtr hostContext = ContextFactory::getContext( Context::Host );
    LAMA_LOG_INFO( logger, "Run local on Host, halo on Host" )
    A.setContext( hostContext, hostContext );
    doVectorTimesMatrixSyncAsyncTests( y, A, x, corResult );
    ContextPtr cudaContext = lama_test::CUDAContext::getContext();
    //2. CUDA, Host
    LAMA_LOG_INFO( logger, "Run local on CUDA, halo on Host" )
    A.setContext( cudaContext, hostContext );
    doVectorTimesMatrixSyncAsyncTests( y, A, x, corResult );
    //3. Host, CUDA
    LAMA_LOG_INFO( logger, "Run local on Host, halo on Cuda" )
    A.setContext( hostContext, cudaContext );
    doVectorTimesMatrixSyncAsyncTests( y, A, x, corResult );
    //4. CUDA, CUDA
    LAMA_LOG_INFO( logger, "Run local on CUDA, halo on CUDA" )
    A.setContext( cudaContext, cudaContext );
    doVectorTimesMatrixSyncAsyncTests( y, A, x, corResult );
    //reset to defaults
    A.setContext( hostContext, hostContext );
}

template<typename MatrixType>
void vectorTimesMatrixTestImpl()
{
    typedef typename MatrixType::MatrixValueType ValueType;
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    CUDAHostContextManager::setAsCurrent( cuda );
    LAMA_LOG_INFO( logger, "set CUDAHostContext as CUDAHostContextManager" )
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    shared_ptr<Distribution> dist( new BlockDistribution( vectorSize, comm ) );
    MatrixType matrixTypeMatrix( dist, dist );
    DenseVector<ValueType> denseVector( dist, 1.0 );
    DenseVector<ValueType> denseTemp( dist, 0.0 );
    // Run vector-matrix multiplication with zero matrix
    doVectorTimesMatrixLocationTests( denseTemp, matrixTypeMatrix, denseVector, denseTemp );
    int numRows = 4 * size;
    int numCols = 4 * size;
    DenseVector<ValueType> denseCorrectResult2( dist, 0.0 );
    LAMAArray<ValueType>& localDenseCorrectResult2 = denseCorrectResult2.getLocalValues();
    scoped_array<ValueType> values( new ValueType[numRows * numCols] );
    {
        HostWriteAccess<ValueType> localDenseCorrectResult2Access( localDenseCorrectResult2 );

        for ( IndexType j = 0; j < numCols; ++j )
        {
            ValueType colSum = 0.0;

            for ( IndexType i = 0; i < numRows; ++i )
            {
                ValueType value = 0.0;

                if ( j == i || j + size == i || j - size == i || j + 2 * size == i || j - 2 * size == i
                        || j + ( numRows - 1 ) == i || j - ( numRows - 1 ) == i )
                {
                    value = 1000.0f * ( i + 1 ) + ( j + 1 );
                }

                values[i * numCols + j] = value;
                colSum += value;
            }

            if ( dist->isLocal( j ) )
            {
                localDenseCorrectResult2Access[dist->global2local( j )] = colSum;
            }
        }
    }
    MatrixType repM;
    repM.setRawDenseData( numRows, numCols, values.get() );
    DenseVector<ValueType> denseVector0( vectorSize, 1.0 );
    DenseVector<ValueType> denseResult0( vectorSize, 0.0 );
    doVectorTimesMatrixLocationTests( denseResult0, repM, denseVector0, denseCorrectResult2 );
    MatrixType matrixTypeMatrix2( repM, dist, dist );
    doVectorTimesMatrixLocationTests( denseTemp, matrixTypeMatrix2, denseVector, denseCorrectResult2 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( vectorTimesMatrixTest, ValueType, test_types )
{
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for CSR" )
    vectorTimesMatrixTestImpl< CSRSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for ELL" )
    vectorTimesMatrixTestImpl< ELLSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for JDS" )
    vectorTimesMatrixTestImpl< JDSSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for COO" )
    vectorTimesMatrixTestImpl< COOSparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for DIA" )
    vectorTimesMatrixTestImpl< DIASparseMatrix<ValueType> >();
    LAMA_LOG_INFO( logger, "vectorTimesMatrixTest: call implementation for Dense" )
    vectorTimesMatrixTestImpl< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */
double randomNumber()
{
    return rand() / static_cast<double>( RAND_MAX );
}

BOOST_AUTO_TEST_CASE( cTorTest )
{
    const IndexType n = 10;
    const float v = 1.0f;
    DenseVector<float> cudaTestVector( n, v );
    const IndexType m = 28576;
    std::vector<float> values( m );
    std::vector<double> valuesD( m );

    for ( IndexType i = 0; i < m; ++i )
    {
        values[i] = static_cast<float>( randomNumber() );
        valuesD[i] = values[i];
    }

    DenseVector<double> hostVector( m, values.data() );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    DenseVector<double> cudaVectorA( m, valuesD.data() );
    DenseVector<double> cudaVectorB( hostVector );
    cudaVectorA.setContext( cuda );
    cudaVectorB.setContext( cuda );

    for ( IndexType i = 0; i < cudaVectorA.size(); i++ )
    {
        BOOST_CHECK_EQUAL( cudaVectorA.getValue( i ), cudaVectorB.getValue( i ) );
    }
}

/* --------------------------------------------------------------------- */

//TODO: Do we need this test? VectorTest/operatorDotProductTest executes the same operations
//with different contexts (e.g. CUDA)
BOOST_AUTO_TEST_CASE( dotProductTest )
{
    IndexType n = 4;
    DenseVector<float> vec1( n, 1.0 );
    DenseVector<float> vec2( n, 2.0 );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    vec1.setContext( cuda );
    vec2.setContext( cuda );
    Scalar result = vec1.dotProduct( vec2 );
    BOOST_CHECK( result == ( 2.0 * n ) );
}

/* --------------------------------------------------------------------- */

//TODO: Do we need this test? VectorTest/SpecialAssignmentTest executes the same operations
//with different contexts (e.g. CUDA)
BOOST_AUTO_TEST_CASE_TEMPLATE( scaleVectorTest, ValueType, test_types )
{
    LAMA_LOG_INFO( logger, "scaleVectorTest<" << common::getScalarType<ValueType>() << ">" )
    IndexType n = 4;
    DenseVector<ValueType> vec1( n, 1.0 );
    DenseVector<ValueType> vec2( n, 0.25 );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    vec1.setContext( cuda );
    vec2.setContext( cuda );
    ValueType alpha = 0.5;
    vec1 *= alpha - 0.25;

    for ( IndexType i = 0; i < vec1.size(); i++ )
    {
        BOOST_CHECK_EQUAL( vec1.getValue( i ), vec2.getValue( i ) );
    }
}

/* --------------------------------------------------------------------- */

//TODO: Do we need this test? VectorTest/AssignmentVectorExpressionTest executes the same operations
//with different contexts (e.g. CUDA)
BOOST_AUTO_TEST_CASE_TEMPLATE( vectorDifferenceTest, ValueType, test_types )
{
    IndexType n = 4;
    DenseVector<ValueType> vec1( n, 3.0 );
    DenseVector<ValueType> vec2( n, 2.0 );
    DenseVector<ValueType> vec3( n, 1.0 );
    DenseVector<ValueType> vec4( n, 0.0 );
    DenseVector<ValueType> vec6( n, 1.0 );
    DenseVector<ValueType> vec7( n, 0.0 );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    vec1.setContext( cuda );
    vec2.setContext( cuda );
    vec3.setContext( cuda );
    vec4.setContext( cuda );
    vec6.setContext( cuda );
    vec7.setContext( cuda );
    MaxNorm maxnorm;
    vec4 = vec1 - vec2;
    DenseVector<ValueType> resultA = vec4 - vec3;
    Scalar resultnormA = maxnorm( resultA );
    BOOST_CHECK( resultnormA == 0.0 );
    DenseVector<ValueType> vec5 = vec1 - vec2;
    DenseVector<ValueType> resultB = vec5 - vec3;
    Scalar resultnormB = maxnorm( resultB );
    BOOST_CHECK( resultnormB == 0.0 );
    vec6 = vec6 - vec3;
    DenseVector<ValueType> resultC = vec6 - vec7;
    Scalar resultnormC = maxnorm( resultC );
    BOOST_CHECK( resultnormC == 0.0 );
}

//TODO: Do we need this test? VectorTest/AssignmentVectorExpressionTest executes the same operations
//with different contexts (e.g. CUDA)

BOOST_AUTO_TEST_CASE_TEMPLATE( scaledVectorsDifferenceTest, ValueType, test_types )
{
    Scalar alpha = 0.5;
    Scalar beta = 0.5;
    MaxNorm maxnorm;
    IndexType n = 4;
    DenseVector<ValueType> vec1( n, 4.0 );
    DenseVector<ValueType> vec2( n, 2.0 );
    DenseVector<ValueType> result( n, 1.0 );
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    vec1.setContext( cuda );
    vec2.setContext( cuda );
    result.setContext( cuda );
    DenseVector<ValueType> vecCalc1( n, 0.0 );
    vecCalc1.setContext( cuda );
    vecCalc1 = alpha * vec1 - beta * vec2;
    DenseVector<ValueType> resultA = vecCalc1 - result;
    Scalar resultnormA = maxnorm( resultA );
    BOOST_CHECK( resultnormA == 0.0 );
    DenseVector<ValueType> vecCalc2( vec1 );
    vecCalc2.setContext( cuda );
    vecCalc2 = alpha * vecCalc2 - beta * vec2;
    DenseVector<ValueType> resultB = vecCalc2 - result;
    Scalar resultnormB = maxnorm( resultB );
    BOOST_CHECK( resultnormB == 0.0 );
    DenseVector<ValueType> vecCalc3( vec2 );
    vecCalc3.setContext( cuda );
    vecCalc3 = alpha * vec1 - beta * vecCalc3;
    DenseVector<ValueType> resultC = vecCalc3 - result;
    Scalar resultnormC = maxnorm( resultC );
    BOOST_CHECK( resultnormC == 0.0 );
    const ValueType gamma = 1.0;
    DenseVector<ValueType> result2( n, 3.0 );
    result2.setContext( cuda );
    DenseVector<ValueType> vecCalc4( n, 6.0 );
    vecCalc4.setContext( cuda );
    vecCalc4 = gamma * vecCalc4 - alpha * vecCalc4;
    DenseVector<ValueType> resultD = vecCalc4 - result2;
    Scalar resultnormD = maxnorm( resultD );
    BOOST_CHECK( resultnormD == 0.0 );
}

/* ------------------------------------------------------------------------- */

template<typename MatrixType>
void operatorMatrixTimeVectorTestMethod()
{
    typedef typename MatrixType::MatrixValueType ValueType;
    MaxNorm maxnorm;
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    {
        IndexType n = 4;
        const ValueType matrixAValues[] =
        { 0.6, 0.0, 0.0, 0.4, 0.7, 0.4, 0.0, 0.0, 0.0, 0.0, 0.9, 0.4, 0.2, 0.5, 0.0, 0.3 };
        CSRSparseMatrix<ValueType> matrixAtmp;
        matrixAtmp.setRawDenseData( n, n, matrixAValues );
        MatrixType matrixA( matrixAtmp );
        matrixA.setContext( cuda );
        DenseVector<ValueType> inputVector( n, 1.0 );
        ValueType resultVectorValues[] =
        { 1.0, 1.1, 1.3, 1.0 };
        DenseVector<ValueType> resultVector( n, resultVectorValues );
        DenseVector<ValueType> calcVector( n, 0.0 );
        calcVector = matrixA * inputVector;
        DenseVector<ValueType> diff = calcVector - resultVector;
        Scalar resultnorm = maxnorm( diff );
        BOOST_CHECK( resultnorm < 1E-5 );
    }
    {
        //Test with an 8x8 Matrix
        IndexType n = 8;
        const CSRSparseMatrix<ValueType> bigMatrixAtmp = TestSparseMatrices::n8m8Laplace1D<ValueType>();
        MatrixType bigMatrixA( bigMatrixAtmp );
        bigMatrixA.setContext( cuda );
        ValueType bigInputVectorValues[] =
        { 1.0, 5.0, 3.0, 2.0, 9.0, 4.0, 7.0, 0.0 };
        DenseVector<ValueType> bigInputVector( n, bigInputVectorValues );
        ValueType bigResultVectorValues[] =
        { -3.0, 6.0, -1.0, -8.0, 12.0, -8.0, 10.0, -7.0 };
        DenseVector<ValueType> bigResultVector( n, bigResultVectorValues );
        DenseVector<ValueType> bigCalcVector( n, 0.0 );
        bigCalcVector = bigMatrixA * bigInputVector;
        DenseVector<ValueType> bigDiff = bigCalcVector - bigResultVector;
        Scalar resultnorm2 = maxnorm( bigDiff );
        BOOST_CHECK( resultnorm2 < 1E-5 );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( operatorSparseMatrixTimeVectorTest, T, test_types )
{
    typedef T ValueType;
    operatorMatrixTimeVectorTestMethod< CSRSparseMatrix<ValueType> >();
    operatorMatrixTimeVectorTestMethod< ELLSparseMatrix<ValueType> >();
    operatorMatrixTimeVectorTestMethod< JDSSparseMatrix<ValueType> >();
    operatorMatrixTimeVectorTestMethod< COOSparseMatrix<ValueType> >();
    operatorMatrixTimeVectorTestMethod< DIASparseMatrix<ValueType> >();
    operatorMatrixTimeVectorTestMethod< DenseMatrix<ValueType> >();
}

/* --------------------------------------------------------------------- */

//TODO: Do we need this test? VectorTest/AssignmentOpMatrixExpressionTest executes the same operations
//with different contexts (e.g. CUDA)
template<typename MatrixType>
void assignmentVectorSubtractExprMatrixTimeVectorTestMethod()
{
    typedef typename MatrixType::MatrixValueType ValueType;
    ContextPtr cuda = lama_test::CUDAContext::getContext();
    {
        EquationHelper::EquationSystem<ValueType> system3x3 = EquationHelper::get3x3SystemA<ValueType>();
        MatrixType matrix3x3( system3x3.coefficients );
        matrix3x3.setContext( cuda );
        DenseVector<ValueType> x3x3( system3x3.solution );
        DenseVector<ValueType> y3x3( system3x3.rhs );
        DenseVector<ValueType> result3x3_0( x3x3.size(), 1.0 );
        result3x3_0 = y3x3 - matrix3x3 * x3x3;
        BOOST_CHECK( maxNorm( result3x3_0 ) < 1e-16 );
        DenseVector<ValueType> result3x3_1 = y3x3 - matrix3x3 * x3x3;
        BOOST_CHECK( maxNorm( result3x3_1 ) < 1e-16 );
        DenseVector<ValueType> result3x3_2 = system3x3.solution;
        result3x3_2 = y3x3 - matrix3x3 * result3x3_2;
        BOOST_CHECK( maxNorm( result3x3_2 ) < 1e-16 );
        DenseVector<ValueType> result3x3_3 = system3x3.rhs;
        result3x3_3 = result3x3_3 - matrix3x3 * x3x3;
        BOOST_CHECK( maxNorm( result3x3_3 ) < 1e-16 );
    }
    {
        EquationHelper::EquationSystem<ValueType> system4x4 = EquationHelper::get4x4SystemA<ValueType>();
        MatrixType matrix4x4( system4x4.coefficients );
        DenseVector<ValueType> x4x4( system4x4.solution );
        DenseVector<ValueType> y4x4( system4x4.rhs );
        DenseVector<ValueType> result4x4_0( x4x4.size(), 1.0 );
        result4x4_0 = y4x4 - matrix4x4 * x4x4;
        BOOST_CHECK( maxNorm( result4x4_0 ) < 1e-16 );
        DenseVector<ValueType> result4x4_1 = y4x4 - matrix4x4 * x4x4;
        BOOST_CHECK( maxNorm( result4x4_1 ) < 1e-16 );
        DenseVector<ValueType> result4x4_2 = system4x4.solution;
        result4x4_2 = y4x4 - matrix4x4 * result4x4_2;
        BOOST_CHECK( maxNorm( result4x4_2 ) < 1e-16 );
        DenseVector<ValueType> result4x4_3 = system4x4.rhs;
        result4x4_3 = result4x4_3 - matrix4x4 * x4x4;
        BOOST_CHECK( maxNorm( result4x4_3 ) < 1e-16 );
        EquationHelper::EquationSystem<ValueType> system8x8 = EquationHelper::get8x8SystemA<ValueType>();
        MatrixType matrix8x8( system8x8.coefficients );
        matrix8x8.setContext( cuda );
        DenseVector<ValueType> x8x8( system8x8.solution );
        DenseVector<ValueType> y8x8( system8x8.rhs );
        DenseVector<ValueType> result8x8_0( x8x8.size(), 1.0 );
        result8x8_0 = y8x8 - matrix8x8 * x8x8;
        BOOST_CHECK( maxNorm( result8x8_0 ) < 1e-16 );
        DenseVector<ValueType> result8x8_1 = y8x8 - matrix8x8 * x8x8;
        BOOST_CHECK( maxNorm( result8x8_1 ) < 1e-16 );
        DenseVector<ValueType> result8x8_2 = system8x8.solution;
        result8x8_2 = y8x8 - matrix8x8 * result8x8_2;
        BOOST_CHECK( maxNorm( result8x8_2 ) < 1e-16 );
        DenseVector<ValueType> result8x8_3 = system8x8.rhs;
        result8x8_3 = result8x8_3 - matrix8x8 * x8x8;
        BOOST_CHECK( maxNorm( result8x8_3 ) < 1e-16 );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignmentVectorSubtractExprSparseMatrixTimeVectorTest, ValueType, test_types )
{
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< CSRSparseMatrix<ValueType> >();
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< ELLSparseMatrix<ValueType> >();
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< DIASparseMatrix<ValueType> >();
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< JDSSparseMatrix<ValueType> >();
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< COOSparseMatrix<ValueType> >();
    assignmentVectorSubtractExprMatrixTimeVectorTestMethod< DenseMatrix<ValueType> >();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
