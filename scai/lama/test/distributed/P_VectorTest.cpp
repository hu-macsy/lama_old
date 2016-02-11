/**
 * @file P_VectorTest.cpp
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
 * @brief Contains the implementation of the class P_VectorTest.
 * @author: Alexander Büchel, Lauretta Schubert
 * @date 06.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <scai/common/unique_ptr.hpp>

#include <scai/lama/Vector.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/dmemo/NoDistribution.hpp>

#include <scai/lama/test/TestSparseMatrices.hpp>
#include <scai/lama/test/EquationHelper.hpp>

#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;
using scai::common::unique_ptr;
using scai::common::scoped_array;
using scai::common::shared_ptr;
using scai::common::Exception;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_VectorTestConfig
{
    P_VectorTestConfig()
    {
        comm = Communicator::getCommunicator();   // default communicator
        m_inputVectorBaseName = scai::test::Configuration::getPath() + "/testVector";
        m_formattedInputVectorBaseName = m_inputVectorBaseName + "Formatted";
        m_xdrDoubleInputVectorBaseName = m_inputVectorBaseName + "XDRDouble";
    }

    ~P_VectorTestConfig()
    {
        comm = CommunicatorPtr();
    }

    std::string m_inputVectorBaseName;
    std::string m_formattedInputVectorBaseName;
    std::string m_xdrDoubleInputVectorBaseName;
};

BOOST_FIXTURE_TEST_SUITE( P_VectorTest, P_VectorTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.P_VectorTest" )

/* --------------------------------------------------------------------- */

//TODO: Decide if this test is neccessary
//BOOST_AUTO_TEST_CASE_TEMPLATE( readWriteTest, ValueType, test_types)
//{
//    if ( 0 != comm->getRank() )
//        return;
//
//    std::string test_formatted_filename = m_formattedInputVectorBaseName+"Tmp";
//    std::string test_binary_filename = m_inputVectorBaseName+"BinaryTmp";
//    std::string test_unformatted_filename = m_inputVectorBaseName
//            +"UnformattedTmp";
//    std::string test_xdr_double_filename = m_inputVectorBaseName+"XDRDoubleTmp";
//    std::string test_xdr_float_filename = m_inputVectorBaseName+"XDTFloatTmp";
//
//    // read ascii and xdr vector files and compare them
//    DenseVector<ValueType> float_ascii_vector(m_formattedInputVectorBaseName);
//    DenseVector<ValueType> float_xdr_vector(m_xdrDoubleInputVectorBaseName);
//    for (IndexType i = 0; i < float_ascii_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( float_xdr_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//
//    // ascii write test
//    float_ascii_vector.writeToFile(test_formatted_filename, sblas::FORMATTED);
//    DenseVector<ValueType> test_float_ascii_vector(test_formatted_filename);
//    for (IndexType i = 0; i < test_float_ascii_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( test_float_ascii_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//
//    // write to binary vector file. read it. compare it with the ascii-file
//    float_ascii_vector.writeToFile(test_binary_filename, sblas::BINARY);
//    DenseVector<ValueType> test_float_binary_vector(test_binary_filename);
//    BOOST_CHECK(==);
//    for (IndexType i = 0; i < test_float_binary_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( test_float_binary_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//
//    // write to binary vector file. read it. compare it with the ascii-file
//    float_ascii_vector.writeToFile(test_unformatted_filename,
//                                   sblas::UNFORMATTED);
//    DenseVector<ValueType>
//            test_float_unformatted_vector(test_unformatted_filename);
//    for (IndexType i = 0; i < test_float_unformatted_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( test_float_unformatted_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//
//    // write a float xdr file
//    float_ascii_vector.writeToFile(test_xdr_float_filename, sblas::XDR,
//                                   sblas::INTERNAL);
//    DenseVector<ValueType> test_float_xdr_vector(test_xdr_float_filename);
//    for (IndexType i = 0; i < test_float_xdr_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( test_float_xdr_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//
//    // write a double xdr file
//    float_ascii_vector.writeToFile(test_xdr_double_filename, sblas::XDR);
//    DenseVector<ValueType> test_double_xdr_vector(test_xdr_double_filename);
//    for (IndexType i = 0; i < float_ascii_vector.size(); ++i )
//    {
//        BOOST_CHECK_EQUAL ( float_ascii_vector.getValue(i) , float_ascii_vector.getValue(i) );
//    }
//}
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void vectorCheck( DenseVector<ValueType>& v, DenseVector<ValueType>& w )
{
    // check equality of two vectors with the same distribution
    IndexType vectorSize = v.size();
    IndexType localVectorSize = v.getLocalValues().size();
    // dense vectors v and w must have the same size and distribution
    BOOST_REQUIRE_EQUAL( vectorSize, w.size() );
    BOOST_REQUIRE_EQUAL( v.getDistribution(), w.getDistribution() );
    BOOST_REQUIRE_EQUAL( localVectorSize, v.getDistribution().getLocalSize() );
    ReadAccess<ValueType> vRead( v.getLocalValues() );
    ReadAccess<ValueType> wRead( w.getLocalValues() );

    // so we just have to compare the local values

    for ( IndexType i = 0; i < localVectorSize; ++i )
    {
        BOOST_CHECK_EQUAL( vRead[i], wRead[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildTest )
{
    typedef DenseVector<double> DoubleVector;
    PartitionId size = comm->getSize();
    IndexType vectorSize = 3 * size;
    scoped_array<double> values( new double[vectorSize] );

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        values[i] = static_cast<double>( i );
    }

    DoubleVector repV( vectorSize, values.get() );

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        BOOST_CHECK_CLOSE( values[i], repV( i ).getValue<double>(), 1e-16 );
    }

    DistributionPtr dist( new BlockDistribution( vectorSize, comm ) );
    DoubleVector distV( repV );
    distV.redistribute( dist );

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        BOOST_CHECK_CLOSE( repV( i ).getValue<double>(), distV( i ).getValue<double>(), 1e-16 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( vectorTimesMatrixTest, ValueType, test_types )
{
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    SCAI_LOG_DEBUG( logger, "VectorSize is " << vectorSize )
    shared_ptr<Distribution> dist( new BlockDistribution( vectorSize, comm ) );
    shared_ptr<Distribution> repdist( new NoDistribution( vectorSize ) );
    CSRSparseMatrix<ValueType> matrixTypeMatrix( dist, repdist );
    matrixTypeMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );
    DenseVector<ValueType> denseVector( dist, 1.0 );
    DenseVector<ValueType> denseResult( repdist, 0.0 );
    const Matrix& matrix = matrixTypeMatrix;
    const Vector& vector = denseVector;
    Vector& result = denseResult;
    SCAI_LOG_INFO( logger, "Vector(NoDist) = Vector(BlockDist) * Matrix(BlockDist, NoDist)" )
    result = vector * matrix;
    ContextPtr host = Context::getHostPtr();
    matrixTypeMatrix.setContextPtr( host, host );

    for ( IndexType i = 0; i < result.size(); ++i )
    {
        BOOST_CHECK_EQUAL( result.getValue( i ), denseResult.getValue( i ) );
    }

    int numRows = 4 * size;
    int numCols = 4 * size;
    DenseVector<ValueType> denseCorrectResult2( dist, 0.0 );
    HArray<ValueType>& localDenseCorrectResult2 =
        denseCorrectResult2.getLocalValues();
    scoped_array<ValueType> values( new ValueType[ numRows * numCols ] );
    {
        WriteAccess<ValueType> localDenseCorrectResult2Access ( localDenseCorrectResult2 );

        for ( IndexType j = 0; j < numCols; ++j )
        {
            ValueType columnSum = 0.0;

            for ( IndexType i = 0; i < numRows; ++i )
            {
                ValueType value = 0.0;

                if ( j == i || j + size == i || j - size == i || j + 2 * size == i || j - 2 * size == i || j + ( numRows - 1 ) == i
                        || j - ( numRows - 1 ) == i )
                {
                    value = static_cast<ValueType>( 1000.0 * ( i + 1 ) + ( j + 1 ) );
                }

                values[ i * numCols + j ] = value;
                columnSum += value;
            }

            if ( dist->isLocal( j ) )
            {
                localDenseCorrectResult2Access[ dist->global2local( j ) ] = columnSum;
            }
        }
    }
    CSRSparseMatrix<ValueType> repM;
    repM.setRawDenseData( numRows, numCols, values.get() );
    repM.setCommunicationKind( Matrix::SYNCHRONOUS );
    DenseVector<ValueType> denseVector0( vectorSize, 1.0 );
    DenseVector<ValueType> denseResult0( vectorSize, 0.0 );
    Vector& result0 = denseResult0;
    SCAI_LOG_INFO( logger, "Vector(rep) = Vector(rep) * Matrix(rep)" )
    result0 = denseVector0 * repM;
    BOOST_CHECK_EQUAL( vectorSize, result.size() );

    for ( IndexType i = 0; i < result.size(); ++i )
    {
        BOOST_CHECK_EQUAL( denseCorrectResult2.getValue( i ), denseResult0.getValue( i ) );
    }

    repM.setContextPtr( host, host );
    CSRSparseMatrix<ValueType> matrixTypeMatrix2( repM, dist, dist );
    const Matrix& matrix2 = matrixTypeMatrix2;
    DenseVector<ValueType> denseVector2( dist, 1.0 );
    DenseVector<ValueType> result2( dist, 0.0 );
    SCAI_LOG_INFO( logger, "Vector(BlockDist) = Vector(BlockDist) * Matrix(BlockDist, BlockDist)" )
    result2 = denseVector2 * matrix2;
    BOOST_CHECK_EQUAL( vectorSize, result2.size() );

    for ( IndexType i = 0; i < result2.size(); ++i )
    {
        BOOST_CHECK_EQUAL( denseCorrectResult2.getValue( i ), result2.getValue( i ) );
    }

    matrixTypeMatrix2.setContextPtr( host, host );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixTimesVectorTest, ValueType, test_types )
{
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    shared_ptr<Distribution> dist( new BlockDistribution( vectorSize, comm ) );
    shared_ptr<Distribution> repdist( new NoDistribution( vectorSize ) );
    CSRSparseMatrix<ValueType> matrixTypeMatrix( dist, repdist );
    DenseVector<ValueType> denseVector( repdist, 1.0 );
    DenseVector<ValueType> denseResult( dist, 0.0 );
    const Matrix& matrix = matrixTypeMatrix;
    const Vector& vector = denseVector;
    Vector& result = denseResult;
    result = matrix * vector;
    ContextPtr host = Context::getHostPtr();
    matrixTypeMatrix.setContextPtr( host, host );

    for ( IndexType i = 0; i < result.size(); ++i )
    {
        BOOST_CHECK_EQUAL( result.getValue( i ), denseResult.getValue( i ) );
    }

    int numRows = 4 * size;
    int numCols = 4 * size;
    DenseVector<ValueType> denseCorrectResult2( dist, 0.0 );
    HArray<ValueType>& localDenseCorrectResult2 =
        denseCorrectResult2.getLocalValues();
    scoped_array<ValueType> values( new ValueType[ numRows * numCols ] );
    {
        WriteAccess<ValueType> localDenseCorrectResult2Access ( localDenseCorrectResult2 );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType rowSum = 0.0;

            for ( IndexType j = 0; j < numCols; ++j )
            {
                ValueType value = 0.0;

                if ( j == i || j + size == i || j - size == i || j + 2 * size == i || j - 2 * size == i || j + ( numRows - 1 ) == i
                        || j - ( numRows - 1 ) == i )
                {
                    value = static_cast<ValueType>( 1000.0 * ( i + 1 ) + ( j + 1 ) );
                }

                values[ i * numCols + j ] = value;
                rowSum += value;
            }

            if ( dist->isLocal( i ) )
            {
                localDenseCorrectResult2Access[ dist->global2local( i ) ] = rowSum;
            }
        }
    }
    CSRSparseMatrix<ValueType> repM;
    repM.setRawDenseData( numRows, numCols, values.get() );
    DenseVector<ValueType> denseVector0( vectorSize, 1.0 );
    DenseVector<ValueType> denseResult0( vectorSize, 0.0 );
    Vector& result0 = denseResult0;
    result0 = repM * denseVector0;

    for ( IndexType i = 0; i < result.size(); ++i )
    {
        BOOST_CHECK_EQUAL( denseCorrectResult2.getValue( i ), denseResult0.getValue( i ) );
    }

    repM.setContextPtr( host, host );
    CSRSparseMatrix<ValueType> matrixTypeMatrix2( repM, dist, dist );
    const Matrix& matrix2 = matrixTypeMatrix2;
    DenseVector<ValueType> denseVector2( dist, 1.0 );
    DenseVector<ValueType> result2( dist, 0.0 );
    result2 = matrix2 * denseVector2;

    for ( IndexType i = 0; i < result2.size(); ++i )
    {
        BOOST_CHECK_EQUAL( denseCorrectResult2.getValue( i ), result2.getValue( i ) );
    }

    matrixTypeMatrix2.setContextPtr( host, host );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignLocalTest, ValueType, test_types )
{
    const IndexType vectorSize = 25;
    shared_ptr<Distribution> dist( new CyclicDistribution( vectorSize, 2, comm ) );
    const IndexType localSize = dist->getLocalSize();
    HArray<float> localData;
// Be careful: for more than 13 processors some of them do not throw exception
    SCAI_CHECK_THROW(
    {   DenseVector<ValueType> denseVector( localData, dist );}, Exception );
    {
        WriteOnlyAccess<float> wLocalData( localData, localSize );

        for ( IndexType i = 0; i < localSize; ++i )
        {
            wLocalData[i] = static_cast<float>( dist->local2global( i ) );
        }
    }
// consider only local constructors / assignments, no redistributions
    DenseVector<ValueType> denseVector1( localData, dist );
    DenseVector<ValueType> denseVector2;
    denseVector2.assign( localData, dist );
    DenseVector<double> denseVector3( denseVector1 );// also local copies
    DenseVector<float> denseVector4;
    denseVector4.assign( denseVector2 );
    BOOST_REQUIRE_EQUAL( denseVector1.size(), vectorSize );
    BOOST_REQUIRE_EQUAL( denseVector2.size(), vectorSize );
    BOOST_REQUIRE_EQUAL( denseVector3.size(), vectorSize );
    BOOST_REQUIRE_EQUAL( denseVector4.size(), vectorSize );
    BOOST_CHECK_EQUAL( denseVector1.getLocalValues().size(), localSize );
    BOOST_CHECK_EQUAL( denseVector2.getLocalValues().size(), localSize );
    BOOST_CHECK_EQUAL( denseVector3.getLocalValues().size(), localSize );
    BOOST_CHECK_EQUAL( denseVector4.getLocalValues().size(), localSize );

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        ValueType expected = static_cast<ValueType>( i );
        Scalar given1 = denseVector1.getValue( i );
        Scalar given2 = denseVector2.getValue( i );
        Scalar given3 = denseVector3.getValue( i );
        Scalar given4 = denseVector4.getValue( i );
        BOOST_CHECK_EQUAL( expected, given1.getValue<ValueType>() );
        BOOST_CHECK_EQUAL( expected, given2.getValue<ValueType>() );
        BOOST_CHECK_EQUAL( expected, given3.getValue<ValueType>() );
        BOOST_CHECK_EQUAL( expected, given4.getValue<ValueType>() );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignValueTest, ValueType, test_types )
{
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    shared_ptr<Distribution> dist( new BlockDistribution( vectorSize, comm ) );
    DenseVector<ValueType> denseVector( dist, 1.0 );
    const Scalar value = 3.0;
    denseVector = value;

    for ( IndexType i = 0; i < denseVector.size(); ++i )
    {
        BOOST_CHECK_EQUAL( value.getValue<ValueType>(), denseVector.getValue( i ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, ValueType, test_types )
{
    const IndexType vectorSize = 100; // global vector size
    const IndexType chunkSize = 3;// used for cyclic distribution
    shared_ptr<Distribution> dist1 ( new BlockDistribution( vectorSize, comm ) );
    shared_ptr<Distribution> dist2 ( new CyclicDistribution( vectorSize, chunkSize, comm ) );
    scoped_array<ValueType> vectorData ( new ValueType[vectorSize] );

    for ( IndexType i = 0; i < vectorSize; i++ )
    {
        vectorData[i] = static_cast<ValueType>( i );
    }

    DenseVector<ValueType> replicatedVector( vectorSize, vectorData.get() );
    SCAI_LOG_INFO( logger, "redistributeTest: distribtedVector ( replicatedVector )" );
// constructor with redistribute that is here a localization
    DenseVector<ValueType> dist1Vector( replicatedVector, dist1 );
    {
        const HArray<ValueType>& localValues = dist1Vector.getLocalValues();
        const IndexType localSize = localValues.size();
        BOOST_REQUIRE_EQUAL( localSize, dist1->getLocalSize() );
        ReadAccess<ValueType> rLocalValues( localValues );

        // each processor checks for correct local values

        for ( IndexType i = 0; i < localSize; ++i )
        {
            ValueType expected = static_cast<ValueType>( dist1->local2global( i ) );
            BOOST_CHECK_EQUAL( expected, rLocalValues[i] );
        }
    }
// redistribute block-distributed vector to cyclic-distributed vector
    SCAI_LOG_INFO( logger, "redistributeTest: distributedVector ( distributedVector, newDist )" );
    DenseVector<double> dist2Vector( dist1Vector );
    dist2Vector.redistribute( dist2 );
    {
        const HArray<double>& localValues = dist2Vector.getLocalValues();
        const IndexType localSize = localValues.size();
        BOOST_REQUIRE_EQUAL( localSize, dist2->getLocalSize() );
        ReadAccess<double> rLocalValues( localValues );

        for ( IndexType i = 0; i < localSize; ++i )
        {
            double expected = dist2->local2global( i );
            BOOST_CHECK_EQUAL( expected, rLocalValues[i] );
        }
    }
    shared_ptr<Distribution> dist3 ( new NoDistribution( vectorSize ) );
    DenseVector<ValueType> dist3Vector( dist2Vector, dist3 );
    dist2Vector.redistribute( dist3 );
    {
        const HArray<double>& localValues2 = dist2Vector.getLocalValues();
        const HArray<ValueType>& localValues3 = dist3Vector.getLocalValues();
        const IndexType localSize = localValues2.size();
        BOOST_REQUIRE_EQUAL( localSize, vectorSize );
        BOOST_REQUIRE_EQUAL( localSize, localValues3.size() );
        ReadAccess<double> rLocalValues2( localValues2 );
        ReadAccess<ValueType> rLocalValues3( localValues3 );

        for ( IndexType i = 0; i < localSize; ++i )
        {
            double expected2 = i;
            ValueType expected3 = static_cast<ValueType>( i );
            BOOST_CHECK_EQUAL( expected2, rLocalValues2[i] );
            BOOST_CHECK_EQUAL( expected3, rLocalValues3[i] );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gatherTest, ValueType, test_types )
{
    PartitionId size = comm->getSize();
    const IndexType vectorSize = 4 * size;
    shared_ptr<Distribution> dist ( new BlockDistribution( vectorSize, comm ) );
    shared_ptr<Distribution> rep ( new NoDistribution( vectorSize ) );
    scoped_array<ValueType> vectorData ( new ValueType[vectorSize] );

    for ( IndexType i = 0; i < vectorSize; i++ )
    {
        vectorData[i] = static_cast<ValueType>( i );
    }

    DenseVector<ValueType> replicatedVector( vectorSize, vectorData.get() );
    DenseVector<ValueType> dist1Vector( replicatedVector, dist );
// gather block-distributed vector to replicated vector
    DenseVector<ValueType> allVector( dist1Vector, rep );
    ReadAccess<ValueType> data( allVector.getLocalValues() );
    BOOST_CHECK_EQUAL ( vectorSize, allVector.getLocalValues().size() );

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        BOOST_CHECK_EQUAL( vectorData[i], data[i] );
    }

    data.release();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ExpressionCtorTest, ValueType, test_types )
{
    // TODO: detect error, when executing with n = 4; for vector times matrix expressions
    PartitionId size = comm->getSize();
    IndexType n = 4 * size;
    CSRSparseMatrix<ValueType> Id (
        TestSparseMatrices::nnIdentityMatrix<ValueType>( n ) );
    DistributionPtr dist( new BlockDistribution( Id.getNumRows(), comm ) );
    Id.redistribute( dist, dist );
    DenseVector<ValueType> testvector5( dist, 5.0 );
    DenseVector<ValueType> testvector3( dist, 3.0 );
    Scalar s2 = 2.0;
    Scalar s3 = 3.0;
    DenseVector<ValueType> resultvectorm1( dist, -1.0 );
    DenseVector<ValueType> resultvector1( dist, 1.0 );
    DenseVector<ValueType> resultvector2( dist, 2.0 );
    DenseVector<ValueType> resultvector7( dist, 7.0 );
    DenseVector<ValueType> resultvector10( dist, 10.0 );
    DenseVector<ValueType> resultvector16( dist, 16.0 );
    DenseVector<ValueType> resultvector13( dist, 13.0 );
    DenseVector<ValueType> resultvector14( dist, 14.0 );
    DenseVector<ValueType> resultvector11( dist, 11.0 );
    DenseVector<ValueType> resultvector19( dist, 19.0 );
    DenseVector<ValueType> resultvector8( dist, 8.0 );
// Vector Expressions in constructors
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector )" );
    DenseVector<ValueType> d3( s2 * testvector5 );
    SCAI_LOG_DEBUG( logger, "dist = " << d3.getDistribution() );
    vectorCheck ( d3, resultvector10 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector + vector )" );
    DenseVector<ValueType> d5( s2 * testvector5 + testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d5.getDistribution() );
    vectorCheck ( d5, resultvector13 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector + scalar * vector)" );
    DenseVector<ValueType> d6( testvector5 + s3 * testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d6.getDistribution() );
    vectorCheck ( d6, resultvector14 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector + scalar * vector)" );
    DenseVector<ValueType> d4( s2 * testvector5 + s3 * testvector3 );
    SCAI_LOG_DEBUG( logger, "d4 = " << d4 );
    vectorCheck ( d4, resultvector19 );
// MatrixVector-Expressions here:
// Note: the constructed vector inherits the row distribution of the matrix
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector )" );
    DenseVector<ValueType> d1( Id * testvector5 );
    SCAI_LOG_DEBUG( logger, "d1 = " << d1 );
    vectorCheck( d1, testvector5 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector )" );
//A*a*x
    DenseVector<ValueType> d17( Id * s2 * testvector5 );
    SCAI_LOG_DEBUG( logger, "dist = " << d17.getDistribution() );
    vectorCheck ( d17, resultvector10 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar )" );
    DenseVector<ValueType> d16( Id * testvector5 * s2 );
    SCAI_LOG_DEBUG( logger, "dist = " << d16.getDistribution() );
    vectorCheck ( d16, resultvector10 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * matrix * vector )" );
    DenseVector<ValueType> d2( s2 * Id * testvector5 );
    SCAI_LOG_DEBUG( logger, "dist = " << d2.getDistribution() );
    vectorCheck ( d2, resultvector10 );
// VectorMatrix-Expressions here:
// Note: the constructed vector inherits the column distribution of the matrix
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix )" );
    DenseVector<ValueType> dd1( testvector5 * Id );
    SCAI_LOG_DEBUG( logger, "d1 = " << d1 );
    vectorCheck( dd1, testvector5 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix)" );
    //A*a*x
    DenseVector<ValueType> dd17( s2 * testvector5 * Id );
    SCAI_LOG_DEBUG( logger, "dist = " << d17.getDistribution() );
    vectorCheck ( dd17, resultvector10 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar )" );
    DenseVector<ValueType> dd16( testvector5 * Id * s2 );
    SCAI_LOG_DEBUG( logger, "dist = " << d16.getDistribution() );
    vectorCheck ( dd16, resultvector10 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix )" );
    DenseVector<ValueType> dd2( s2 * testvector5 * Id );
    SCAI_LOG_DEBUG( logger, "dist = " << d2.getDistribution() );
    vectorCheck ( dd2, resultvector10 );
// Plus-MatrixVectorExpressions:
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector + vector )" );
    DenseVector<ValueType> d11( Id * testvector5 + testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d11.getDistribution() );
    vectorCheck ( d11, resultvector8 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector + scalar * vector )" );
    DenseVector<ValueType> d7( Id * testvector5 + s2 * testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d7.getDistribution() );
    vectorCheck ( d7, resultvector11 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector + vector * scalar )" );
    DenseVector<ValueType> d9( Id * testvector5 + testvector3 * s2 );
    vectorCheck ( d9, resultvector11 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector + scalar * vector )" );
    DenseVector<ValueType> d8( Id * s2 * testvector5 + s2 * testvector3 );
    vectorCheck ( d8, resultvector16 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar + scalar * vector )" );
    DenseVector<ValueType> d10( Id * testvector5 * s2 + s2 * testvector3 );
    vectorCheck ( d10, resultvector16 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector + vector )" );
    DenseVector<ValueType> d13( Id * s2 * testvector5 + testvector3 );
    vectorCheck ( d13, resultvector13 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * matrix * vector + vector )" );
    DenseVector<ValueType> d14( s2 * Id * testvector5 + testvector3 );
    vectorCheck ( d14, resultvector13 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * matrix * vector + scalar * vector )" );
    DenseVector<ValueType> d15( s2 * Id * testvector5 + s3 * testvector3 );
    vectorCheck ( d15, resultvector19 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar + vector )" );
    DenseVector<ValueType> d12( Id * testvector5 * s2 + testvector3 );
    vectorCheck ( d12, resultvector13 );
// VectorMatrix Plus Vector Expressions:
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix + vector )" );
    DenseVector<ValueType> d111( testvector5 * Id + testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d11.getDistribution() );
    vectorCheck ( d111, resultvector8 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix + scalar * vector )" );
    DenseVector<ValueType> d71( testvector5 * Id + s2 * testvector3 );
    SCAI_LOG_DEBUG( logger, "dist = " << d7.getDistribution() );
    vectorCheck ( d71, resultvector11 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix + vector * scalar )" );
    DenseVector<ValueType> d91( testvector5 * Id + testvector3 * s2 );
    vectorCheck ( d91, resultvector11 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix + scalar * vector )" );
    DenseVector<ValueType> d81( s2 * testvector5 * Id + s2 * testvector3 );
    vectorCheck ( d81, resultvector16 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar + scalar * vector )" );
    DenseVector<ValueType> d101( testvector5 * Id * s2 + s2 * testvector3 );
    vectorCheck ( d101, resultvector16 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector  * matrix+ vector )" );
    DenseVector<ValueType> d131( s2 * testvector5 * Id + testvector3 );
    vectorCheck ( d131, resultvector13 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix + vector )" );
    DenseVector<ValueType> d141( s2 * testvector5 * Id + testvector3 );
    vectorCheck ( d141, resultvector13 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix + scalar * vector )" );
    DenseVector<ValueType> d151( s2 * testvector5 * Id + s3 * testvector3 );
    vectorCheck ( d151, resultvector19 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar + vector )" );
    DenseVector<ValueType> d121( testvector5 * Id * s2 + testvector3 );
    vectorCheck ( d121, resultvector13 );
//Minus Vector Expressions:
    SCAI_LOG_INFO( logger, "Ctor: vector( vector - vector )" );
    DenseVector<ValueType> d18( testvector5 - testvector3 );
    vectorCheck ( d18, resultvector2 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector - vector )" );
    DenseVector<ValueType> d19( s2 * testvector5 - testvector3 );
    vectorCheck ( d19, resultvector7 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector - scalar * vector )" );
    DenseVector<ValueType> d20( testvector5 - s2 * testvector3 );
    vectorCheck ( d20, resultvectorm1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector - scalar * vector )" );
    DenseVector<ValueType> d21( s2 * testvector5 - s3 * testvector3 );
    vectorCheck ( d21, resultvector1 );
//Minus MatrixVector Expressions:
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector - scalar * vector )" );
    DenseVector<ValueType> d22( Id * s2 * testvector5 - s3 * testvector3 );
    vectorCheck ( d22, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar - scalar * vector )" );
    DenseVector<ValueType> d23( Id * testvector5 * s2 - s3 * testvector3 );
    vectorCheck ( d23, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector - vector * scalar )" );
    DenseVector<ValueType> d24( Id * s2 * testvector5 - testvector3 * s3 );
    vectorCheck ( d24, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar - vector * scalar )" );
    DenseVector<ValueType> d25( Id * testvector5 * s2 - testvector3 * s3 );
    vectorCheck ( d25, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector - vector )" );
    DenseVector<ValueType> d26( Id * testvector5 - testvector3 );
    vectorCheck ( d26, resultvector2 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * vector * scalar - scalar * vector )" );
    DenseVector<ValueType> d27( Id * testvector5 * s2 - testvector3 );
    vectorCheck ( d27, resultvector7 );
    SCAI_LOG_INFO( logger, "Ctor: vector( matrix * scalar * vector - vector )" );
    DenseVector<ValueType> d28( Id * s2 * testvector5 - testvector3 );
    vectorCheck ( d28, resultvector7 );
//VectorMatrix Minus Vector Expressions:
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix - scalar * vector )" );
    DenseVector<ValueType> d222( s2 * testvector5 * Id - s3 * testvector3 );
    vectorCheck ( d222, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar - scalar * vector )" );
    DenseVector<ValueType> d232( testvector5 * Id * s2 - s3 * testvector3 );
    vectorCheck ( d232, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix - vector * scalar )" );
    DenseVector<ValueType> d242( s2 * testvector5 * Id - testvector3 * s3 );
    vectorCheck ( d242, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar - vector * scalar )" );
    DenseVector<ValueType> d252( testvector5 * Id * s2 - testvector3 * s3 );
    vectorCheck ( d252, resultvector1 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix - vector )" );
    DenseVector<ValueType> d262( testvector5 * Id - testvector3 );
    vectorCheck ( d262, resultvector2 );
    SCAI_LOG_INFO( logger, "Ctor: vector( vector * matrix * scalar - scalar * vector )" );
    DenseVector<ValueType> d272( testvector5 * Id * s2 - testvector3 );
    vectorCheck ( d272, resultvector7 );
    SCAI_LOG_INFO( logger, "Ctor: vector( scalar * vector * matrix - vector )" );
    DenseVector<ValueType> d282( s2 * testvector5 * Id - testvector3 );
    vectorCheck ( d282, resultvector7 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ExpressionAssignmentOperatorTest, ValueType, test_types )
{
    // TODO: detect error, when executing with n = 4; for vector times matrix expressions
    PartitionId size = comm->getSize();
    IndexType n = 4 * size;
    CSRSparseMatrix<ValueType> Id (
        TestSparseMatrices::nnIdentityMatrix<ValueType>( n ) );
    DistributionPtr dist( new BlockDistribution( Id.getNumRows(), comm ) );
    Id.redistribute( dist, dist );
    DenseVector<ValueType> vector( dist, 1.0 );
    DenseVector<ValueType> testvector5( dist, 5.0 );
    DenseVector<ValueType> testvector3( dist, 3.0 );
    Scalar s2 = 2.0;
    Scalar s3 = 3.0;
    DenseVector<ValueType> resultvectorm1( dist, -1.0 );
    DenseVector<ValueType> resultvector1( dist, 1.0 );
    DenseVector<ValueType> resultvector2( dist, 2.0 );
    DenseVector<ValueType> resultvector7( dist, 7.0 );
    DenseVector<ValueType> resultvector8( dist, 8.0 );
    DenseVector<ValueType> resultvector10( dist, 10.0 );
    DenseVector<ValueType> resultvector11( dist, 11.0 );
    DenseVector<ValueType> resultvector13( dist, 13.0 );
    DenseVector<ValueType> resultvector14( dist, 14.0 );
    DenseVector<ValueType> resultvector16( dist, 16.0 );
    DenseVector<ValueType> resultvector19( dist, 19.0 );
// make sure that vectorCheck works correctly
    vectorCheck( resultvector1, resultvector1 );
//VectorExpressions:
// vector = scalar * vector1
    vector = s2 * testvector5;
    vectorCheck ( vector, resultvector10 );
//a*x+y
    vector = s2 * testvector5 + testvector3;
    vectorCheck ( vector, resultvector13 );
//x+s*y
    vector = testvector5 + s3 * testvector3;
    vectorCheck ( vector, resultvector14 );
//a*x+b*y
    vector = s2 * testvector5 + s3 * testvector3;
    vectorCheck ( vector, resultvector19 );
//MatrixVector-Expressions:
//A*x
    vector = Id * testvector5;
    vectorCheck( vector, testvector5 );
//A*a*x
    vector = Id * s2 * testvector5;
    vectorCheck ( vector, resultvector10 );
//A*x*a
    vector = Id * testvector5 * s2;
    vectorCheck ( vector, resultvector10 );
//a*A*x
    vector = s2 * Id * testvector5;
    vectorCheck ( vector, resultvector10 );
//VectorMatrix-Expressions:
    //x*A
    vector = testvector5 * Id;
    vectorCheck( vector, testvector5 );
    //a*x*A
    vector = s2 * testvector5 * Id;
    vectorCheck ( vector, resultvector10 );
    //x*a*A
    vector = testvector5 * s2 * Id;
    vectorCheck ( vector, resultvector10 );
    //x*A*a
    vector = testvector5 * Id * s2;
    vectorCheck ( vector, resultvector10 );
//MatrixVectorPlusVector Expressions:
//A*x+y
    vector = Id * testvector5 + testvector3;
    vectorCheck ( vector, resultvector8 );
//A*x+a*y
    vector = Id * testvector5 + s2 * testvector3;
    vectorCheck ( vector, resultvector11 );
//A*x+y*a
    vector = Id * testvector5 + testvector3 * s2;
    vectorCheck ( vector, resultvector11 );
//A*a*x+b*y
    vector = Id * s2 * testvector5 + s2 * testvector3;
    vectorCheck ( vector, resultvector16 );
//A*x*a+b*y
    vector = Id * testvector5 * s2 + s2 * testvector3;
    vectorCheck ( vector, resultvector16 );
//A*a*x+y
    vector = Id * s2 * testvector5 + testvector3;
    vectorCheck ( vector, resultvector13 );
//a*A*x+y
    vector = s2 * Id * testvector5 + testvector3;
    vectorCheck ( vector, resultvector13 );
//a*A*x+b*y
    vector = s2 * Id * testvector5 + s3 * testvector3;
    vectorCheck ( vector, resultvector19 );
//A*x*a+b*y
    vector = Id * testvector5 * s2 + testvector3;
    vectorCheck ( vector, resultvector13 );
//VectorMatrixPlusVector Expressions:
    //x*A+y
    vector = testvector5 * Id + testvector3;
    vectorCheck ( vector, resultvector8 );
    //x*A+a*y
    vector = testvector5 * Id + s2 * testvector3;
    vectorCheck ( vector, resultvector11 );
    //x*A+y*a
    vector = testvector5 * Id + testvector3 * s2;
    vectorCheck ( vector, resultvector11 );
    //a*x*A+b*y
    vector = s2 * testvector5 * Id + s2 * testvector3;
    vectorCheck ( vector, resultvector16 );
    //x*a*A+b*y
    vector = testvector5 * s2 * Id + s2 * testvector3;
    vectorCheck ( vector, resultvector16 );
    //x*A*a+y
    vector = testvector5 * Id * s2 + testvector3;
    vectorCheck ( vector, resultvector13 );
    //a*x*A+y
    vector = s2 * testvector5 * Id + testvector3;
    vectorCheck ( vector, resultvector13 );
    //x*a*A+b*y
    vector = testvector5 * s2 * Id + s3 * testvector3;
    vectorCheck ( vector, resultvector19 );
    //x*A*a+b*y
    vector = testvector5 * Id * s2 + testvector3;
    vectorCheck ( vector, resultvector13 );
//Minus Vector Expressions:
//x-y
    vector = testvector5 - testvector3;
    vectorCheck ( vector, resultvector2 );
//a*x-y
    vector = s2 * testvector5 - testvector3;
    vectorCheck ( vector, resultvector7 );
//x-b*y
    vector = testvector5 - s2 * testvector3;
    vectorCheck ( vector, resultvectorm1 );
//a*x-b*y
    vector = s2 * testvector5 - s3 * testvector3;
    vectorCheck ( vector, resultvector1 );
//MatrixVector Minus Vector Expressions:
//A*a*x-b*y
    vector = Id * s2 * testvector5 - s3 * testvector3;
    vectorCheck ( vector, resultvector1 );
//A*x*a-b*y
    vector = Id * testvector5 * s2 - s3 * testvector3;
    vectorCheck ( vector, resultvector1 );
//A*x*a-y*b
    vector = Id * s2 * testvector5 - testvector3 * s3;
    vectorCheck ( vector, resultvector1 );
//A*x*a-y*b
    vector = Id * testvector5 * s2 - testvector3 * s3;
    vectorCheck ( vector, resultvector1 );
//A*x-y
    vector = Id * testvector5 - testvector3;
    vectorCheck ( vector, resultvector2 );
//A*x*a-y*b
    vector = Id * testvector5 * s2 - testvector3;
    vectorCheck ( vector, resultvector7 );
//A*a*x-y*b
    vector = Id * s2 * testvector5 - testvector3;
    vectorCheck ( vector, resultvector7 );
//VectorMatrix Minus Vector Expressions:
    //a*x*A-b*y
    vector = s2 * testvector5 * Id - s3 * testvector3;
    vectorCheck ( vector, resultvector1 );
    //x*a*A-b*y
    vector = testvector5 * s2 * Id - s3 * testvector3;
    vectorCheck ( vector, resultvector1 );
    //x*A*a-y*b
    vector = s2 * Id * testvector5 - testvector3 * s3;
    vectorCheck ( vector, resultvector1 );
    //a*x*A-y*b
    vector = s2 * testvector5 * Id - testvector3 * s3;
    vectorCheck ( vector, resultvector1 );
    //x*A-y
    vector = testvector5 * Id - testvector3;
    vectorCheck ( vector, resultvector2 );
    //x*A*a-y*b
    vector = testvector5 * Id * s2 - testvector3;
    vectorCheck ( vector, resultvector7 );
    //a*A*x-y*b
    vector = s2 * Id * testvector5 - testvector3;
    vectorCheck ( vector, resultvector7 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( dotProductTest, ValueType, test_types )
{
    const IndexType n = 8;
    DistributionPtr dist( new BlockDistribution( n, comm ) );
    DistributionPtr dist1( new BlockDistribution( 8, comm ) );
    DistributionPtr dist2( new CyclicDistribution( 8, 3, comm ) );
    DenseVector<ValueType> v1( dist, 1.0 );
    DenseVector<ValueType> v2( dist, 2.0 );
    DenseVector<ValueType> v3( dist, -2.0 );
    DenseVector<ValueType> v4( dist2, 14.0 );
    Scalar result;

    //Should throw exception, because of different distributions.
    //vectors are distributed if np > 1
    if ( comm->getSize() > 1 )
    {
        SCAI_CHECK_THROW( { result = v3.dotProduct( v4 ); }, Exception );
    }
    else
    {
        SCAI_LOG_INFO( logger, "dotProductTest did not run, because np = 1" );
    }
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
