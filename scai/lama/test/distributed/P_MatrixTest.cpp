/**
 * @file P_MatrixTest.cpp
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
 * @brief Tests for all kind of distributed matrices
 * @author: Thomas Brandes
 * @date 26.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/distribution/CyclicDistribution.hpp>
#include <scai/lama/distribution/GeneralDistribution.hpp>
#include <scai/lama/distribution/GenBlockDistribution.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/lama/test/TestSparseMatrices.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/test/SameMatrixHelper.hpp>
#include <scai/lama/test/TestSparseMatrices.hpp>

using namespace scai;
using namespace scai::lama;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

DistributionPtr makeDistribution( const IndexType n, CommunicatorPtr comm, int kind )
{
    if ( kind == 0 )
    {
        return DistributionPtr( new BlockDistribution( n, comm ) );
    }
    else if ( kind == 1 )
    {
        return DistributionPtr( new CyclicDistribution( n, 3, comm ) );
    }
    else if ( kind == 2 )
    {
        float weight = 1.0f;

        if ( comm->getRank() % 2 == 1 )
        {
            weight = 0.0001f;
        }

        return DistributionPtr( new GenBlockDistribution( n, weight, comm ) );
    }
    else if ( kind == 3 )
    {
        float weight = 1.0f;

        if ( comm->getRank() % 2 == 0 )
        {
            weight = 0.0001f;
        }

        return DistributionPtr( new GenBlockDistribution( n, weight, comm ) );
    }
    else if ( kind == 4 )
    {
        PartitionId size = comm->getSize();
        PartitionId rank = comm->getRank();
        IndexType val = 1713;
        std::vector<IndexType> localIndexes;

        for ( IndexType i = 0; i < n; i++ )
        {
            PartitionId owner = val % size;

            if ( owner == rank )
            {
                localIndexes.push_back( i );
            }

            val = ( val * 21 ) % 1913;
        }

        return DistributionPtr( new GeneralDistribution( n, localIndexes, comm ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported kind of makeDistribution" );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( P_MatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixTest" );

typedef boost::mpl::list <
CSRSparseMatrix<float>,
                DenseMatrix<double>,
                DIASparseMatrix<float>,
                ELLSparseMatrix<double>  > MatrixTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( LocalConstructorTest, MatrixType, MatrixTypes )
{
    CommunicatorPtr comm = Communicator::getCommunicator();
    const int size = comm->getSize();
    const int n = 10;
    DistributionPtr dist( new BlockDistribution( n * size, comm ) );
    CSRStorage<float> csrGlobal; // global replicated storage
    csrGlobal.setIdentity( n * size );
    BOOST_REQUIRE_EQUAL( csrGlobal.getNumRows(), n * size );
    CSRStorage<double> csrLocal;// local distributed storage
    csrLocal.localize( csrGlobal, *dist );
    BOOST_REQUIRE_EQUAL( csrLocal.getNumRows(), n );
    MatrixType m1;
    m1.assign( csrGlobal );// full replicated matrix
    MatrixType m2 ( csrGlobal, dist, dist );
    MatrixType m3 ( csrLocal, dist, dist );
    BOOST_CHECK_EQUAL( m2.getDistribution(), m3.getDistribution() );
    BOOST_CHECK_EQUAL( m2.getColDistribution(), m3.getColDistribution() );
    testSameMatrix( m1, m2 );
    testSameMatrix( m1, m3 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetDenseDataTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    CommunicatorPtr comm = Communicator::getCommunicator();
    // use multiple test cases for different nValues
    const IndexType nValues[] = { 0, 1, 3, 17 };
    const int nCases = sizeof( nValues ) / sizeof( IndexType );

    for ( int k = 0; k < nCases; ++k )
    {
        int n = nValues[k];
        DistributionPtr dist( new BlockDistribution( n, comm ) );
        MatrixType m1;
        m1.setIdentity( dist );
        // build CSR data owned by this partition for identity matrix
        std::vector<ValueType> values;
        values.reserve( dist->getLocalSize() * n );
        MatrixType m2;   // will be filled with dense data
        ValueType eps = 0.05;

        // eps will only be used for sparse matrices

        if ( m2.getMatrixKind() == Matrix::DENSE )
        {
            eps = 0.0;
        }

        for ( IndexType i = 0; i < n; ++i )
        {
            if ( ! dist->isLocal( i ) )
            {
                continue;
            }

            for ( IndexType j = 0; j < n; ++j )
            {
                ValueType val = static_cast<ValueType>( eps );

                if ( i == j )
                {
                    val = static_cast<ValueType>( 1 );
                }

                values.push_back( val );
            }
        }

        m2.setRawDenseData( dist, dist, &values[0], 2 * eps );
        testSameMatrix( m1, m2 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetCSRDataTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    CommunicatorPtr comm = Communicator::getCommunicator();
    // use multiple test cases for different nValues
    const IndexType nValues[] = { 0, 1, 3, 17 };
    const int nCases = sizeof( nValues ) / sizeof( IndexType );

    for ( int k = 0; k < nCases; ++k )
    {
        int n = nValues[k];
        DistributionPtr dist( new BlockDistribution( n, comm ) );
        MatrixType m1;
        m1.setIdentity( dist );
        // build CSR data owned by this partition for identity matrix
        std::vector<IndexType> iaOffset;
        std::vector<IndexType> iaSizes;
        std::vector<IndexType> ja;
        std::vector<ValueType> values;
        IndexType offset = 0;
        IndexType numValues = 0;
        iaOffset.push_back( offset );

        for ( IndexType i = 0; i < n; ++i )
        {
            if ( ! dist->isLocal( i ) )
            {
                continue;
            }

            ja.push_back( i );
            ValueType val = static_cast<ValueType>( 1.0 );
            values.push_back( val );
            offset ++;
            numValues ++;
            iaOffset.push_back( offset );
            iaSizes.push_back( 1 );
        }

        MatrixType m2;
        m2.setRawCSRData( dist, dist, numValues, &iaOffset[0], &ja[0], &values[0] );
        testSameMatrix( m1, m2 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CopyConstructorTest, MatrixType, MatrixTypes )
{
    const IndexType N1 = 10;
    const IndexType N2 = 10;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    CSRSparseMatrix<double> inputA;
    MatrixCreator<double>::buildPoisson2D( inputA, 9, N1, N2 );
    CommunicatorPtr comm = Communicator::getCommunicator();
    const IndexType n = inputA.getNumRows();
    SCAI_LOG_DEBUG( logger, "inputA = " << inputA )
    DistributionPtr dist( new BlockDistribution( n, comm ) );
    MatrixType m1( inputA );
    SCAI_LOG_DEBUG( logger, "m1( inputA ) = " << m1 )
    testSameMatrix( inputA, m1 );
    MatrixType m2( inputA, dist, dist );
    SCAI_LOG_DEBUG( logger, "m2( inputA, dist, dist ) = " << m2 )
    testSameMatrix( inputA, m2 );
    MatrixType m3;
    m3 = inputA;
    SCAI_LOG_DEBUG( logger, "m3 ( = inputA ) = " << m3 )
    testSameMatrix( inputA, m3 );
};

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list<CSRSparseMatrix<float>, DIASparseMatrix<double>, JDSSparseMatrix<double>, COOSparseMatrix<float>,
        ELLSparseMatrix<double> > SparseMatrixTypes;

BOOST_AUTO_TEST_CASE_TEMPLATE( SwapTest, MatrixType, SparseMatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType globalSize = 100;
    StorageType localStorage;
    CommunicatorPtr comm = Communicator::getCommunicator();
    DistributionPtr dist( new BlockDistribution( globalSize, comm ) );
    DistributionPtr rep ( new NoDistribution( globalSize ) );
    const IndexType localSize = dist->getLocalSize();
    localStorage.setIdentity( localSize );
    MatrixType unityMatrix( dist, dist ); // not yet unity
    unityMatrix.swapLocalStorage( localStorage );// now unity matrix
    MatrixType nullMatrix( dist, dist );
    localStorage.setIdentity( localSize );
    localStorage.setDiagonal( 0.0 );// null matrix with diagonal property
    nullMatrix.swapLocalStorage( localStorage );
    nullMatrix.swap( unityMatrix );// now we have swapped
// get the diagonals
    DenseVector<ValueType> nullVector;
    DenseVector<ValueType> oneVector;
    nullMatrix.getDiagonal( oneVector );
    unityMatrix.getDiagonal( nullVector );
// now compare
    nullVector.redistribute( rep );
    oneVector.redistribute( rep );
    BOOST_REQUIRE_EQUAL( nullVector.getLocalValues().size(), globalSize );
    BOOST_REQUIRE_EQUAL( oneVector.getLocalValues().size(), globalSize );
    ReadAccess<ValueType> rNull( nullVector.getLocalValues() );
    ReadAccess<ValueType> rOne( oneVector.getLocalValues() );
    ValueType null = 0.0;
    ValueType one = 1.0;

    for ( IndexType i = 0; i < globalSize; ++i )
    {
        BOOST_CHECK_EQUAL( rNull[i], null );
        BOOST_CHECK_EQUAL( rOne[i], one );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FullConstructorTest, MatrixType, SparseMatrixTypes )
{
    std::string prefix = scai::test::Configuration::getPath();
    SCAI_LOG_INFO( logger, "prefix = " << prefix << ", can be changed by LAMA" );
    CSRSparseMatrix<double> tmp( prefix + "/can___24.mtx" );
    SCAI_LOG_INFO( logger, "constructed replicated matrix by file: " << tmp );
    IndexType numTotalRows = tmp.getNumRows();
    CommunicatorPtr comm = Communicator::getCommunicator();
    DistributionPtr dist = makeDistribution( numTotalRows, comm, 1 );
    SCAI_LOG_INFO( logger, "new distribution: " << *dist );
    CSRSparseMatrix<double> matrix( tmp, dist, dist );
    IndexType numLocalRows = matrix.getLocalNumRows();

    /* get distributed data */

    IndexType numLocalValues = matrix.getLocalStorage().getNumValues();
    IndexType numHaloValues = matrix.getHaloStorage().getNumValues();

    SCAI_LOG_INFO( logger, *comm << ": local N = " << numLocalRows << ", nnz = "
                   << numLocalValues << " (local) + " << numHaloValues << " (halo)" );

    common::scoped_array<IndexType> iaLocal( new IndexType[ numLocalRows + 1 ] );
    common::scoped_array<IndexType> jaLocal( new IndexType[ numLocalValues ] );
    common::scoped_array<double> valuesLocal( new double[ numLocalValues ] );

    common::scoped_array<IndexType> iaHalo( new IndexType[ numLocalRows + 1 ] );
    common::scoped_array<IndexType> jaHalo( new IndexType[ numHaloValues ] );
    common::scoped_array<double> valuesHalo( new double[ numHaloValues ] );

    const CSRStorage<double>& localSt = matrix.getLocalStorage();

    ReadAccess<IndexType> iaLocalRead( localSt.getIA() );
    ReadAccess<IndexType> jaLocalRead( localSt.getJA() );
    ReadAccess<double> valuesLocalRead( localSt.getValues() );
    const CSRStorage<double>& haloSt = matrix.getHaloStorage();
    ReadAccess<IndexType> iaHaloRead( haloSt.getIA() );
    ReadAccess<IndexType> jaHaloRead( haloSt.getJA() );
    ReadAccess<double> valuesHaloRead( haloSt.getValues() );
    ReadAccess<IndexType> halo2global( matrix.getHalo().getRequiredIndexes() );

    for ( IndexType i = 0; i < numLocalRows + 1; ++i )
    {
        iaLocal[i] = iaLocalRead[i];
        SCAI_LOG_TRACE( logger, "local: ia[ " << i << " ] = " << iaLocal[i] )
    }

    // Be careful, halo might be 0 x 0, so we have now iaHalo

    if ( haloSt.getNumRows() == 0 )
    {
        SCAI_LOG_TRACE( logger, "local: ia[ 0 .. " << numLocalRows << " + 1 ] = 0 " )

        for ( IndexType i = 0; i < numLocalRows + 1; ++i )
        {
            iaHalo[i] = 0;
        }
    }
    else
    {
        BOOST_REQUIRE_EQUAL( numLocalRows, haloSt.getNumRows() );

        for ( IndexType i = 0; i < numLocalRows + 1; ++i )
        {
            iaHalo[i] = iaHaloRead[i];
            SCAI_LOG_TRACE( logger, "halo: ia[ " << i << " ] = " << iaHalo[i] )
        }
    }

    for ( IndexType i = 0; i < numLocalValues; ++i )
    {
        jaLocal[i] = jaLocalRead[i];   // local indexes remain local
        valuesLocal[i] = valuesLocalRead[i];
        SCAI_LOG_TRACE( logger, "local: ja[ " << i << " ] = " << jaLocal[i] )
        SCAI_LOG_TRACE( logger, "local: values[ " << i << " ] = " << valuesLocal[i] )
    }

    for ( IndexType i = 0; i < numHaloValues; ++i )
    {
        jaHalo[i] = jaHaloRead[i];
        jaHalo[i] = halo2global[jaHalo[i]];   // halo indexes must be global
        valuesHalo[i] = valuesHaloRead[i];
        SCAI_LOG_TRACE( logger, "halo: ja[ " << i << " ] = " << jaHalo[i] )
        SCAI_LOG_TRACE( logger, "halo: values[ " << i << " ] = " << valuesHalo[i] )
    }

    std::vector<IndexType> globalIndexes( numLocalRows );

    for ( IndexType i = 0; i < numLocalRows; ++i )
    {
        globalIndexes[i] = dist->local2global( i );
        SCAI_LOG_TRACE( logger, *comm << ": local row " << i << " is global row " << globalIndexes[i] );
    }

    MatrixType distMatrix( numLocalRows, numLocalValues, numHaloValues, 
                           iaLocal.get(), jaLocal.get(), valuesLocal.get(), 
                           iaHalo.get(), jaHalo.get(), valuesHalo.get(), 
                           globalIndexes, comm );

    testSameMatrix( matrix, distMatrix );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( InvertTest, MatrixType, SparseMatrixTypes )
{
    MatrixType m1 = TestSparseMatrices::n4m4TestMatrix1<double>();
    CommunicatorPtr comm = Communicator::getCommunicator();
    const IndexType n = m1.getNumRows();
    DistributionPtr bdist( new BlockDistribution( n, comm ) );
    DistributionPtr ndist( new NoDistribution( n ) );
    m1.redistribute( bdist, ndist );
    SCAI_LOG_TRACE( logger, "Input random matrix for invert: " << m1 );
    MatrixType m2;
    m2.invert( m1 );
    SCAI_LOG_INFO( logger, "Inverted matrix: " << m2 );
    SCAI_ASSERT_EQUAL_ERROR( m2.getDistribution(), m1.getDistribution() );
    SCAI_ASSERT_EQUAL_ERROR( m2.getColDistribution(), m1.getColDistribution() );
    m1.redistribute( bdist, bdist );// otherwise no matrix mult possible
    // dense matrix multiplication not supported yet with distributed m2
    MatrixType mm( m1 * m2 );
    SCAI_LOG_INFO( logger, "Result of matrix x inverse : " << mm );
    MatrixType unity;
    unity.setIdentity( bdist );
    SCAI_LOG_INFO( logger, "Distributed identity matrix: " << mm );
    // CLOSE test not sufficient as mm might have inexact ZERO values, so use small value
    Scalar small( common::TypeTraits<typename MatrixType::MatrixValueType>::small() );
    testSameMatrix( unity, mm, small );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixMultTest, MatrixType, SparseMatrixTypes )
{
    typedef typename MatrixType::MatrixValueType ValueType;
// first we test for replicated matrices
    MatrixType matrix1 = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    MatrixType matrix2 = TestSparseMatrices::n4m3MatrixE2<ValueType>();
    MatrixType matrixP( matrix1 * matrix2 );
    MatrixType matrixR = TestSparseMatrices::n6m3MatrixERes<ValueType>();
    SCAI_LOG_INFO( logger, "verify: n6m4MatrixE1 * n4m3MatrixE2 = n6m3MatrixDRes" );
    testSameMatrix( matrixR, matrixP );
// now we use distributions
    CommunicatorPtr comm = Communicator::getCommunicator();
    DistributionPtr rowDist( new BlockDistribution( matrix1.getNumRows(), comm ) );
    DistributionPtr colDist( new BlockDistribution( matrix1.getNumColumns(), comm ) );
    DistributionPtr repDist( new NoDistribution( matrix2.getNumColumns() ) );
    matrix1.redistribute ( rowDist, colDist );
    matrix2.redistribute ( colDist, repDist );
    matrixP = matrix1 * matrix2;
    BOOST_REQUIRE_EQUAL( matrixP.getDistribution(), *rowDist );
    BOOST_REQUIRE_EQUAL( matrixP.getColDistribution(), *repDist );
    testSameMatrix( matrixR, matrixP );
    matrixR.redistribute( rowDist, repDist );
// P and R have same distribution, use maxDiffNorm
    Scalar diff = matrixR.maxDiffNorm( matrixP );
    BOOST_CHECK_SMALL( diff.getValue<double>() , 1e-5 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

