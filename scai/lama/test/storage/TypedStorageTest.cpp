/**
 * @file TypedStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Test cases applied to each typed storage class, i.e. test (virtual) methods of MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/logging.hpp>

#include <memory>

using namespace scai;
using namespace lama;
using common::Exception;

using hmemo::HArray;
using utilskernel::HArrayUtils;

using boost::test_tools::per_element;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( TypedStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.TypedStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( factoryTest, ValueType, scai_numeric_test_types )
{
    TypedStorages<ValueType> allMatrixStorages;    // is created by factory
    size_t nFormats = static_cast<size_t>( Format::UNDEFINED );
    nFormats -= 1;  // stencil not in Factory
    SCAI_LOG_INFO( logger, "factoryTest<" << common::TypeTraits<ValueType>::id() << "> : "
                   << allMatrixStorages.size() << " storages"                        )
    BOOST_CHECK_EQUAL( nFormats, allMatrixStorages.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( purgeTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "purgeTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );   // in any case not a non-zero
        SCAI_LOG_DEBUG( logger, "purgeTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        size_t usedBytes = storage.getMemoryUsage();
        storage.clear();   // resets sizes to zero but does not delete any data
        size_t usedBytesCleared = storage.getMemoryUsage();
        storage.purge();   // purge will free allocated data
        size_t usedBytesPurged = storage.getMemoryUsage();
        SCAI_LOG_DEBUG( logger, "purgeTest, purged storage = " << storage << ", used bytes = " << usedBytes
                        << ", used bytes (clear) = " << usedBytesCleared
                        << ", used bytes (purge) = " << usedBytesPurged )
        BOOST_CHECK( usedBytesCleared <= usedBytes );
        BOOST_CHECK( usedBytesPurged < usedBytes );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( normTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "normTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    DenseStorage<ValueType> dense;
    setDenseData( dense );
    const HArray<ValueType> denseData = dense.getData();

    auto expectedL1Norm = HArrayUtils::l1Norm( denseData );
    auto expectedL2Norm = HArrayUtils::l2Norm( denseData );
    auto expectedMaxNorm = HArrayUtils::maxNorm( denseData );

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_DEBUG( logger, "normTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        setDenseData( storage );
        RealType<ValueType> maxNorm = storage.maxNorm();
        BOOST_CHECK_CLOSE( expectedMaxNorm, maxNorm, 1 );
        RealType<ValueType> l1Norm = storage.l1Norm();
        BOOST_CHECK_CLOSE( expectedL1Norm, l1Norm, 1 );
        RealType<ValueType> l2Norm = storage.l2Norm();
        BOOST_CHECK_CLOSE( expectedL2Norm, l2Norm, 1 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_numeric_test_types )
{
    DenseStorage<ValueType> cmpStorage;
    setDenseData( cmpStorage );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );
        SCAI_LOG_DEBUG( logger, "scaleTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        storage.scale( 2.0 );  // should be executed on context

        for ( IndexType i = 0; i < storage.getNumRows(); ++i )
        {
            for ( IndexType j = 0; j < storage.getNumColumns(); ++j )
            {
                BOOST_CHECK_EQUAL( 2.0 * cmpStorage.getValue( i, j ), storage.getValue( i, j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conjTest, ValueType, scai_numeric_test_types )
{
    // Use random data for the matrix to have complex numbers

    const IndexType numRows    = 3;
    const IndexType numColumns = 3;

    // Initialize dense array, with sparse random values

    HArray<ValueType> denseValues( numRows * numColumns, ValueType( 0) );
    HArrayUtils::fillRandom( denseValues, 1, 0.2f );

    HArray<ValueType> x( numColumns );  // initialization not necessray

    HArrayUtils::fillRandom( x, 1, 1.0f ); // completely random values between 0 and 1

    HArray<ValueType> xconj;
    HArrayUtils::unaryOp( xconj, x, common::UnaryOp::CONJ );

    const ValueType alpha = 1.0;
    const ValueType beta  = 0.0;
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        // y1 = storage * x;
        // y2 = ( storage.conj * x.conj ). conj
        // Proof: y1 == y2
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_INFO( logger, "conj test for storage : " << storage << " on " << *context );
        storage.assign( DenseStorage<ValueType>( numRows, numColumns, denseValues ) );
        HArray<ValueType> z( storage.getNumRows(), ValueType( 0 ), context );
        HArray<ValueType> y1( storage.getNumRows(), ValueType( 0 ), context );
        HArray<ValueType> y2( storage.getNumRows(), ValueType( 0 ), context );
        storage.matrixTimesVector( y1, alpha, x, beta, z, common::MatrixOp::NORMAL );
        storage.conj();
        storage.matrixTimesVector( y2, alpha, xconj, beta, z, common::MatrixOp::NORMAL );
        HArrayUtils::unaryOp( y2, y2, common::UnaryOp::CONJ );
        // y1 == y2 ( per element ) is not valid here
        BOOST_CHECK( HArrayUtils::maxDiffNorm( y1, y2 ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( diagonalTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseSquareData( storage );  // should be square
        SCAI_LOG_DEBUG( logger, "diagonalTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        HArray<ValueType> diag( context );
        storage.getDiagonal( diag );
        BOOST_CHECK_EQUAL( diag.size(), storage.getNumRows() ); // square matrix
        SCAI_LOG_INFO( logger, "diagonalTest: get diagonal = " << diag )
        {
            hmemo::ReadAccess<ValueType> readDiag ( diag );

            for ( IndexType i = 0; i < diag.size(); ++i )
            {
                ValueType s = storage.getValue( i, i );
                BOOST_CHECK_EQUAL( s, readDiag[i] );
            }
        }
        storage.setDiagonal( 0 );

        for ( IndexType i = 0; i < diag.size(); ++i )
        {
            ValueType s = storage.getValue( i, i );
            BOOST_CHECK_EQUAL( s, 0 );
        }

        storage.setDiagonalV( diag );
        {
            hmemo::ReadAccess<ValueType> readDiag ( diag );

            for ( IndexType i = 0; i < diag.size(); ++i )
            {
                ValueType s = storage.getValue( i, i );
                BOOST_CHECK_EQUAL( s, readDiag[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( assignDiagonalTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    hmemo::HArray<ValueType> diagonal( { 3, 2, 1, 4 } );

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
 
        setDenseData( storage );   // do some intialization just to make sure that values are reset

        storage.assignDiagonal( diagonal );

        BOOST_REQUIRE_EQUAL( storage.getNumRows(), diagonal.size() );
        BOOST_REQUIRE_EQUAL( storage.getNumColumns(), diagonal.size() );

        auto csr = convert<CSRStorage<ValueType>>( storage );

        BOOST_REQUIRE_EQUAL( csr.getNumValues(), diagonal.size() );

        BOOST_TEST( hostReadAccess( diagonal ) == hostReadAccess( csr.getValues() ),
                    boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setRowTest, ValueType, scai_numeric_test_types )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    TypedStorages<ValueType> allMatrixStorages( context );

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages assign dense data" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        setDenseData( storage );

        HArray<ValueType> row;

        for ( IndexType i = 0; i < storage.getNumRows(); ++i )
        {
            storage.getRow( row, i );
            storage.setRow( row, i, common::BinaryOp::SUB );
        }

        BOOST_CHECK( storage.maxNorm() < RealType( 0.0001 ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setColumnTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    TypedStorages<ValueType> allMatrixStorages( context );

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages assign dense data" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        setDenseData( storage );

        HArray<ValueType> column;

        for ( IndexType j = 0; j < storage.getNumColumns(); ++j )
        {
            storage.getColumn( column, j );
            storage.setColumn( column, j, common::BinaryOp::SUB );
        }

        // as get/set column is using same data type, storage must now be completely zero

        BOOST_CHECK_EQUAL( storage.maxNorm(), 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( inverseTestIdentity )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    const IndexType N = 10;
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        storage.setIdentity( N );
        std::unique_ptr<MatrixStorage<ValueType> > inverse( storage.newMatrixStorage() );
        inverse->invert( storage );
        BOOST_CHECK( inverse->maxDiffNorm( storage ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( inverseTestRandom )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseRandom( storage );
        BOOST_REQUIRE_EQUAL( storage.getNumRows(), storage.getNumColumns () );  // must be square
        std::unique_ptr<MatrixStorage<ValueType> > inverse( storage.newMatrixStorage() );
        inverse->invert( storage );
        std::unique_ptr<MatrixStorage<ValueType> > result( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > identity( storage.newMatrixStorage() );
        // result = storage * inverse, must be identity
        result->matrixTimesMatrix( ValueType( 1 ), storage, *inverse, ValueType( 0 ), *inverse );
        identity->setIdentity( storage.getNumRows() );
        SCAI_LOG_DEBUG( logger, "max diff norm = " << result->maxDiffNorm( *identity ) )
        BOOST_CHECK( result->maxDiffNorm( *identity ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( inverseTestRandom1 )
{
    // same as inverseTestRandom, but here we invert in place
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseRandom( storage );
        std::unique_ptr<MatrixStorage<ValueType> > inverse( storage.copy() );
        inverse->invert( *inverse );
        std::unique_ptr<MatrixStorage<ValueType> > result( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > identity( storage.newMatrixStorage() );
        result->matrixTimesMatrix( ValueType( 1 ), storage, *inverse, ValueType( 0 ), *inverse );
        identity->setIdentity( storage.getNumRows() );
        SCAI_LOG_DEBUG( logger, "max diff norm = " << result->maxDiffNorm( *identity ) )
        BOOST_CHECK( result->maxDiffNorm( *identity ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixTimesVectorTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    // test of result = alpha * storage * x + beta * y
    // case 1: general
    // case 2: beta = 0, y is zero array
    // case 3: beta = 1, y == result, result += alpha * storage * x + beta
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    const ValueType alpha = 1;
    const ValueType beta  = 2;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 0.5;
    DenseStorage<ValueType> denseStorage( context );
    setRandomData( denseStorage, 3, 4 );
    const HArray<ValueType> x( denseStorage.getNumColumns(), xVal );
    const HArray<ValueType> y( denseStorage.getNumRows(), yVal );
    const HArray<ValueType> yDummy( 0, yVal );
    // for comparison: denseResult = alpha * DenseStorage * x + beta * y
    HArray<ValueType> denseResult1;
    HArray<ValueType> denseResult2;
    HArray<ValueType> denseResult3( y );
    HArray<ValueType> denseResult4;
    denseStorage.matrixTimesVector( denseResult1, alpha, x, beta, y, common::MatrixOp::NORMAL );
    denseStorage.matrixTimesVector( denseResult2, alpha, x, 0, yDummy, common::MatrixOp::NORMAL );
    denseStorage.matrixTimesVector( denseResult3, alpha, x, 1, denseResult3, common::MatrixOp::NORMAL );
    utilskernel::HArrayUtils::compute( denseResult4, beta, common::BinaryOp::MULT, y, context );
    SCAI_LOG_INFO( logger, "matrixTimesVectorTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        storage.assign( denseStorage );
        SCAI_LOG_DEBUG( logger, "GEMV: storage = " << storage )
        // result = alpha * storage * x + beta * y
        HArray<ValueType> result1( context );
        HArray<ValueType> result2( context );
        HArray<ValueType> result3;
        HArrayUtils::assign( result3, y, context );
        // wrong sized y, should throw exception in any case
        BOOST_CHECK_THROW(
        {
            storage.matrixTimesVector( result1, alpha, x, beta, yDummy, common::MatrixOp::NORMAL );
        }, common::AssertException );
        storage.matrixTimesVector( result1, alpha, x, beta, y, common::MatrixOp::NORMAL );
        storage.matrixTimesVector( result2, alpha, x, 0, yDummy, common::MatrixOp::NORMAL );
        storage.matrixTimesVector( result3, alpha, x, 1, result3, common::MatrixOp::NORMAL );
        // should be the same as computed with dense storage
        ValueType eps = 0.001;
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult1, result1 ) < eps );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult2, result2 ) < eps );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult3, result3 ) < eps );
        HArray<ValueType> result4( context );
        storage.allocate( y.size(), x.size() );  // numRows x numCols, all zero
        // test multiplication with zero matrix
        storage.matrixTimesVector( result4, alpha, x, beta, y, common::MatrixOp::NORMAL );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult4, result4 ) < eps );
        // result5 = 0 * storage * x + beta * y, storage can be anything as not used at all
        HArray<ValueType> result5( context );
        storage.clear();  // with alpha == 0, storage is not used at all
        storage.matrixTimesVector( result5, ValueType( 0 ), x, beta, y, common::MatrixOp::NORMAL );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult4, result5 ) < eps );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixTimesVectorAsyncTest )
{
    // test of result = alpha * storage * x + beta * y, asynchronous, general case only
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    const ValueType alpha = 1;
    const ValueType beta  = 2;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 0.5;
    DenseStorage<ValueType> denseStorage;
    setDenseData( denseStorage );
    const HArray<ValueType> x( denseStorage.getNumColumns(), xVal );
    const HArray<ValueType> y( denseStorage.getNumRows(), yVal );
    // for comparison: denseResult = alpha * DenseStorage * x + beta * y
    HArray<ValueType> denseResult;
    denseStorage.matrixTimesVector( denseResult, alpha, x, beta, y, common::MatrixOp::NORMAL );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "matrixTimesVectorTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );
        SCAI_LOG_DEBUG( logger, "storage = " << storage )
        // result = alpha * storage * x + beta * y
        HArray<ValueType> result( context );
        std::unique_ptr<tasking::SyncToken> token( storage.matrixTimesVectorAsync( result, alpha, x, beta, y, common::MatrixOp::NORMAL ) );
        token->wait();
        // should be the same as computed with dense storage
        BOOST_CHECK_EQUAL( HArrayUtils::maxDiffNorm( denseResult, result ), 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( vectorTimesMatrixTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    // test of result = alpha * storage * x + beta * y
    // case 1: general
    // case 2: beta = 0, y is zero array
    // case 3: beta = 1, y == result, result += alpha * storage * x + beta
    const ValueType alpha = 1;
    const ValueType beta  = 2;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 0.5;
    DenseStorage<ValueType> denseStorage;
    setRandomData( denseStorage, 5, 8 );
    const HArray<ValueType> x( denseStorage.getNumRows(), xVal );
    const HArray<ValueType> y( denseStorage.getNumColumns(), yVal );
    const HArray<ValueType> yDummy( 0, yVal );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    // for comparison: denseResult = alpha * x * DenseStorage + beta * y
    HArray<ValueType> denseResult1;
    HArray<ValueType> denseResult2;
    HArray<ValueType> denseResult3( y );
    denseStorage.matrixTimesVector( denseResult1, alpha, x, beta, y, common::MatrixOp::TRANSPOSE );
    denseStorage.matrixTimesVector( denseResult2, alpha, x, 0, yDummy, common::MatrixOp::TRANSPOSE );
    denseStorage.matrixTimesVector( denseResult3, alpha, x, 1, denseResult3, common::MatrixOp::TRANSPOSE );
    HArray<ValueType> denseResult4;
    HArrayUtils::compute( denseResult4, y, common::BinaryOp::MULT, beta );  // result4 = beta * y, for zero matrix
    SCAI_LOG_INFO( logger, "matrixTransposedTimesVector<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        storage.assign( denseStorage );
        SCAI_LOG_DEBUG( logger, "storage = " << storage )
        HArray<ValueType> result1( context );
        HArray<ValueType> result2( context );
        HArray<ValueType> result3;
        HArrayUtils::assign( result3, y, context );
        HArray<ValueType> result4( context );
        storage.matrixTimesVector( result1, alpha, x, beta, y, common::MatrixOp::TRANSPOSE );
        storage.matrixTimesVector( result2, alpha, x, 0, yDummy, common::MatrixOp::TRANSPOSE );
        storage.matrixTimesVector( result3, alpha, x, 1, result3, common::MatrixOp::TRANSPOSE );
        ValueType eps = 0.001;
        // should be the same as computed with dense storage
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult1, result1 ) < eps );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult2, result2 ) < eps );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult3, result3 ) < eps );
        storage.allocate( x.size(), y.size() );  // numRows x numCols, all zero
        // test multiplication with zero matrix
        storage.matrixTimesVector( result4, alpha, x, beta, y, common::MatrixOp::TRANSPOSE );
        BOOST_CHECK( HArrayUtils::maxDiffNorm( denseResult4, result4 ) < eps );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixTimesVectorSparseTest )
{
    // test of y += alpha * storage, where only some rows are filled in storage
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    const ValueType alpha = 2;
    const ValueType beta  = 1;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 2.1;
    DenseStorage<ValueType> denseStorage;
    setDenseHalo( denseStorage );
    const HArray<ValueType> x( denseStorage.getNumColumns(), xVal );
    HArray<ValueType> denseY( denseStorage.getNumRows(), yVal );
    // for comparison: denseY += alpha * DenseStorage * x
    denseStorage.matrixTimesVector( denseY, alpha, x, beta, denseY, common::MatrixOp::NORMAL );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "matrixTimesVectoriSparseTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_DEBUG( logger, "GEMV sparse storage = " << storage << ", set dense halo data" )
        storage.setCompressThreshold( 1.0f );  // introduce row indexes in any case
        setDenseHalo( storage );
        Format format = storage.getFormat();

        if ( format == Format::CSR || format == Format::ELL )
        {
            // these storage format should have sparse row indexes
            BOOST_CHECK( storage.getRowIndexes().size() > 0 );
        }

        HArray<ValueType> y( denseStorage.getNumRows(), yVal, context );
        // y += alpha * storage * x , where only some rows are filled
        storage.matrixTimesVector( y, alpha, x, beta, y, common::MatrixOp::NORMAL );

        // should be the same as computed with dense storage
        BOOST_TEST( hostReadAccess( denseY ) == hostReadAccess( y ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( vectorTimesMatrixAsyncTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    // test of result = alpha * storage * x + beta * y, asynchronous, general case
    // case 1: general
    // case 2: beta = 0, y is zero array
    // case 3: beta = 1, y == result, result += alpha * storage * x + beta
    const ValueType alpha = 1;
    const ValueType beta  = 2;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 0.5;
    DenseStorage<ValueType> denseStorage;
    setDenseData( denseStorage );
    const HArray<ValueType> x( denseStorage.getNumRows(), xVal );
    const HArray<ValueType> y( denseStorage.getNumColumns(), yVal );
    // for comparison: denseResult = alpha * x * DenseStorage + beta * y
    HArray<ValueType> denseResult;
    denseStorage.matrixTimesVector( denseResult, alpha, x, beta, y, common::MatrixOp::TRANSPOSE );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "vectorTimesMatrixAsyncTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );
        SCAI_LOG_DEBUG( logger, "GEVM asynchron: storage = " << storage )
        HArray<ValueType> result( context );
        std::unique_ptr<tasking::SyncToken> token;
        token.reset( storage.matrixTimesVectorAsync( result, alpha, x, beta, y, common::MatrixOp::TRANSPOSE ) );
        token->wait();
        // should be the same as computed with dense storage
        BOOST_CHECK_EQUAL( HArrayUtils::maxDiffNorm( denseResult, result ), 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixVectorTimesSparseTest )
{
    // test of y += alpha * x * storage, where only some rows are filled in storage
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    const ValueType alpha = 2;
    const ValueType beta  = 1;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 2.1;
    DenseStorage<ValueType> denseStorage;
    setDenseHalo( denseStorage );
    const HArray<ValueType> x( denseStorage.getNumRows(), xVal );
    HArray<ValueType> denseY( denseStorage.getNumColumns(), yVal );
    // for comparison: denseY += alpha * x * DenseStorage
    denseStorage.matrixTimesVector( denseY, alpha, x, beta, denseY, common::MatrixOp::TRANSPOSE );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "matrixVecotrTimesSparseTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_DEBUG( logger, "GEMV sparse storage = " << storage << ", set dense halo data" )
        storage.setCompressThreshold( 1.0f );  // introduce row indexes in any case
        setDenseHalo( storage );
        Format format = storage.getFormat();

        if ( format == Format::CSR || format == Format::ELL )
        {
            // these storage format should have sparse row indexes
            BOOST_CHECK( storage.getRowIndexes().size() > 0 );
        }

        HArray<ValueType> y( denseStorage.getNumColumns(), yVal, context );
        // y += alpha * storage * x , where only some rows are filled
        storage.matrixTimesVector( y, alpha, x, beta, y, common::MatrixOp::TRANSPOSE );
        // should be the same as computed with dense storage

        BOOST_TEST( hostReadAccess( denseY ) == hostReadAccess( y ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( transposeTest )
{
    // storage1 = transpose( storage2 ), verify storage1 * x = x * storage 2, for all formats each storage
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "transposeTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages1( context );    // storage for each storage format
    TypedStorages<ValueType> allMatrixStorages2( context );    // storage for each storage format

    for ( size_t s1 = 0; s1 < allMatrixStorages1.size(); ++s1 )
    {
        MatrixStorage<ValueType>& storage1 = *allMatrixStorages1[s1];
        setDenseData( storage1 );
        SCAI_LOG_DEBUG( logger, "storage1 " << s1 << " of " << allMatrixStorages1.size() << " =  " << storage1 )
        // we do a matrix multiplication y = A * x with this storage
        HArray<ValueType> x( context );
        HArray<ValueType> y( storage1.getNumRows(), ValueType( 0 ) );
        HArray<ValueType> result1( context );
        x.resize( storage1.getNumColumns() );
        utilskernel::HArrayUtils::setRandom( x, 1 );
        ValueType alpha = 1;
        ValueType beta  = 0;
        storage1.matrixTimesVector( result1, alpha, x, beta, result1, common::MatrixOp::NORMAL );

        for ( size_t s2 = 0; s2 < allMatrixStorages1.size(); ++s2 )
        {
            MatrixStorage<ValueType>& storage2 = *allMatrixStorages2[s2];
            SCAI_LOG_DEBUG( logger, "storage1 " << s1 << " of " << allMatrixStorages1.size() << " =  " << storage1
                            << ", assign transpose to storage2 " << s2 << " = " << storage2 );
            storage2.assignTranspose( storage1 );
            BOOST_CHECK_EQUAL( storage1.getNumRows(), storage2.getNumColumns() );
            BOOST_CHECK_EQUAL( storage2.getNumRows(), storage1.getNumColumns() );
            HArray<ValueType> result2( context );
            HArray<ValueType> y(  storage2.getNumRows(), ValueType( 0 ) );
            storage2.matrixTimesVector( result2, alpha, x, beta, result2, common::MatrixOp::TRANSPOSE );
            // results should be the same
            BOOST_CHECK( HArrayUtils::maxDiffNorm( result1, result2 ) < common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( jacobiTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "jacobiTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseSquareData( storage );
        // must be square matrix
        SCAI_ASSERT_EQUAL( storage.getNumRows(), storage.getNumColumns(), "storage not square" );
        SCAI_LOG_DEBUG( logger, "storage for jacobiIterate = " << storage )
        const HArray<ValueType> oldSolution( storage.getNumRows(), 1 );
        const HArray<ValueType> rhs( storage.getNumRows(), 2 );
        // clone the storage and set its diagonal to zero, but keep inverse of diagonal
        std::unique_ptr<MatrixStorage<ValueType> > storage1( storage.copy() );
        HArray<ValueType> diagonalInverse;
        storage1->getDiagonal( diagonalInverse );
        HArrayUtils::compute( diagonalInverse, ValueType( 1 ), common::BinaryOp::DIVIDE, diagonalInverse );
        storage1->setDiagonal( 0 );

        ValueType omegas[] = { 1.0, 0.8, 0.5 };
        const int NCASES = sizeof( omegas ) / sizeof( ValueType );

        for ( int k = 0; k < NCASES; ++k )
        {
            ValueType omega = omegas[k];

            SCAI_LOG_DEBUG( logger, "run jacobi, omega = " << omega << ", format = " << storage.getFormat() 
                                    << ", type = " << storage.getValueType() )

            HArray<ValueType> solution1( context );
            HArray<ValueType> solution2( context );
            storage.jacobiIterate( solution1, oldSolution, rhs, omega );
            const ValueType alpha = -1;
            const ValueType beta  = 1;
            //  solution2 = omega * ( rhs - B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
            storage1->matrixTimesVector( solution2, alpha, oldSolution, beta, rhs, common::MatrixOp::NORMAL );
            HArrayUtils::binaryOp( solution2, solution2, diagonalInverse, common::BinaryOp::MULT );
            utilskernel::HArrayUtils::arrayPlusArray( solution2, omega, solution2,
                    ValueType( 1 ) - omega, oldSolution, solution2.getValidContext() );
            // solution1 and solution2 must be the same

            BOOST_CHECK( HArrayUtils::maxDiffNorm( solution1, solution2 ) < common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( jacobiAsyncTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "jacobiTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseSquareData( storage );
        // must be square matrix
        SCAI_ASSERT_EQUAL( storage.getNumRows(), storage.getNumColumns(), "storage not square" );
        SCAI_LOG_DEBUG( logger, "storage for jacobiIterate = " << storage )
        const HArray<ValueType> oldSolution( storage.getNumRows(), 1 );
        const HArray<ValueType> rhs( storage.getNumRows(), 2 );
        ValueType omegas[] = { 1.0, 0.8, 0.5 };
        const int NCASES = sizeof( omegas ) / sizeof( ValueType );

        for ( int k = 0; k < NCASES; ++k )
        {
            ValueType omega = omegas[k];

            SCAI_LOG_DEBUG( logger, "run jacobi async, omega = " << omega << ", format = " << storage.getFormat() 
                                    << ", type = " << storage.getValueType() )

            HArray<ValueType> solution1( context );
            HArray<ValueType> solution2( context );

            storage.jacobiIterate( solution1, oldSolution, rhs, omega );
            {
                std::unique_ptr<tasking::SyncToken> token;
                token.reset( storage.jacobiIterateAsync( solution2, oldSolution, rhs, omega ) );
            }
            BOOST_CHECK( HArrayUtils::maxDiffNorm( solution1, solution2 ) < common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( jacobiHaloTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "jacobiHaloTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        if ( storage.getFormat() == Format::DIA )
        {
            continue;   // DIA  has still bug for diagonal property
        }

        setDenseHalo( storage );
        std::unique_ptr<MatrixStorage<ValueType> > local( storage.newMatrixStorage() );
        setDenseSquareData( *local );
        SCAI_LOG_DEBUG( logger, "storage for jacobiIterateHalo = " << storage )
        HArray<ValueType> diagonal;
        const HArray<ValueType> oldSolution( storage.getNumColumns(), 1 );
        // clone the storage and set its diagonal to zero, but keep inverse of diagonal
        std::unique_ptr<MatrixStorage<ValueType> > storage1( storage.copy() );
        HArray<ValueType> diagonalInverse;
        local->getDiagonal( diagonal );
        HArrayUtils::compute( diagonalInverse, ValueType( 1 ), common::BinaryOp::DIVIDE, diagonal );
        storage1->scaleRows( diagonalInverse );
        ValueType omegas[] = { 1.0, 0.8, 0.5 };
        const int NCASES = sizeof( omegas ) / sizeof( ValueType );

        for ( int k = 0; k < NCASES; ++k )
        {
            ValueType omega = omegas[k];
            HArray<ValueType> solution1( storage.getNumRows(), 1, context );
            HArray<ValueType> solution2( storage.getNumRows(), 1, context );
            // solution1 -= omega * ( B(halo) * oldSolution ) * dinv
            storage.jacobiIterateHalo( solution1, diagonal, oldSolution, omega );
            const ValueType alpha = -omega;
            const ValueType beta  = 1;
            storage1->matrixTimesVector( solution2, alpha, oldSolution, beta , solution2, common::MatrixOp::NORMAL );
            // now solution1 and solution2 must be the same
            /* for debug
            std::cout << "Solution1 :";
            for ( IndexType i = 0; i < storage.getNumRows(); ++i )
            {
                std::cout << " ";
                std::cout << solution1[i];
            }
            std::cout << std::endl;

            std::cout << "Solution2 :";
            for ( IndexType i = 0; i < storage.getNumRows(); ++i )
            {
                std::cout << " ";
                std::cout << solution2[i];
            }
            std::cout << std::endl;
            std::cout << "max diff norm = " << solution1.maxDiffNorm( solution2 );
            std::cout << ", small = " << common::TypeTraits<ValueType>::small() << std::endl;
            */

            BOOST_CHECK( HArrayUtils::maxDiffNorm( solution1, solution2 ) < common::TypeTraits<ValueType>::small() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixAddTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    const IndexType N = 3;
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "matrixAddTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_DEBUG( logger, "storage for matrixAdd = " << storage )
        std::unique_ptr<MatrixStorage<ValueType> > a( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > b( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > c( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > res( storage.newMatrixStorage() );
        ValueType alpha = 1;
        ValueType beta  = 2;
        a->setIdentity( N );
        b->setIdentity( N );
        c->setIdentity( N );
        res->setIdentity( N );
        a->setDiagonal( 2 );
        b->setDiagonal( 3 );
        c->setDiagonal( -3 );
        storage.matrixTimesMatrix( alpha, *a, *b, beta, *c );
        res->setDiagonal( alpha * 2 * 3 - beta * 3 );
        BOOST_REQUIRE_EQUAL( N, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumColumns() );
        BOOST_CHECK_EQUAL( 0, storage.maxDiffNorm( *res ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( binaryOpTest )
{
    // only one ValueType as we here just test correct handling of sparse formats

    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    //   Example data
    //
    //    1  0  1   2       2  1  0   2        -1  -1  1 0
    //    2  1  3   0   -   3  0  -1  2    =   -1   1  4  -2
    //    1  0  2   0       1  1  0   0         0   -1  2  0

    const IndexType m = 3;
    const IndexType n = 4;

    HArray<ValueType> data1     ( {  1,  0, 1, 2,  2, 1,  3,  0, 1,  0, 2, 0 } );
    HArray<ValueType> data2     ( {  2,  1, 0, 2,  3, 0, -1,  2, 1,  1, 0, 0 } );
    HArray<ValueType> resultSub ( { -1, -1, 1, 0, -1, 1,  4, -2, 0, -1, 2, 0 } );
    HArray<ValueType> resultAdd ( {  3,  1, 1, 4,  5, 1,  2,  2, 2,  1, 2, 0 } );

    DenseStorage<ValueType> storage1( m, n, data1 );
    DenseStorage<ValueType> storage2( m, n, data2 );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t i = 0; i < allMatrixStorages.size(); ++i )
    {
        auto& mStorage2 = *allMatrixStorages[i];

        // Case 1: with dense storage

        mStorage2.assign( storage2 );

        mStorage2.binaryOp( storage1, common::BinaryOp::SUB, mStorage2 );

        auto mStorageR = convert<DenseStorage<ValueType>>( mStorage2 );

        BOOST_TEST( hostReadAccess( resultSub ) == hostReadAccess( mStorageR.getValues() ), per_element() );

        // Case 2: with sprse storage

        mStorage2.assign( storage2 );

        auto sparseStorage1 = convert<CSRStorage<ValueType>>( storage1 );

        mStorage2.binaryOp( sparseStorage1, common::BinaryOp::ADD, mStorage2 );

        mStorageR.assign( mStorage2 );

        BOOST_TEST( hostReadAccess( resultAdd ) == hostReadAccess( mStorageR.getValues() ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixMultTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    SCAI_LOG_INFO( logger, "matrixMultTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_DEBUG( logger, "storage for matrixMult = " << storage )
        std::unique_ptr<MatrixStorage<ValueType> > a( storage.newMatrixStorage() );
        std::unique_ptr<MatrixStorage<ValueType> > b( storage.newMatrixStorage() );
        ValueType alpha = 1;
        ValueType beta  = 0;
        setDenseSquareData( *a );
        setDenseSquareData( *b );
        storage.matrixTimesMatrix( alpha, *a, *b, beta, storage );
        HArray<ValueType> x( a->getNumColumns(), 1 );
        HArray<ValueType> dummy;
        HArray<ValueType> y1;
        HArray<ValueType> y2;
        HArray<ValueType> tmp;
        // compute y1 = a * b * x and y2 = storage * x, y1 and y2 must be same
        b->matrixTimesVector( tmp, alpha, x, beta, dummy, common::MatrixOp::NORMAL );
        a->matrixTimesVector( y1, alpha, tmp, beta, dummy, common::MatrixOp::NORMAL );
        storage.matrixTimesVector( y2, alpha, x, beta, dummy, common::MatrixOp::NORMAL );

        BOOST_CHECK( HArrayUtils::maxDiffNorm( y1, y2 ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( symmetryTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );
        SCAI_LOG_DEBUG( logger, "symmetryTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        BOOST_CHECK( !storage.checkSymmetry() );
        setSymDenseData( storage );
        BOOST_CHECK( storage.checkSymmetry() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setCSRDataTest )
{
    using namespace hmemo;
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    ContextPtr context = Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        std::unique_ptr<MatrixStorage<ValueType> > storageDense( storage.copy() );
        SCAI_LOG_DEBUG( logger, "setCSRData for " << storage )
        storage.clear();
        IndexType numRows;
        IndexType numColumns;
        HArray<IndexType> matrixRowSizes( context );
        HArray<IndexType> matrixJA( context );
        HArray<ValueType> matrixValues( context );
        HArray<ValueType> matrixDense( context );
        getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
        // Now we can use matrixRowSizes as IA array
        storage.setCSRData( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues );
        SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                       << ", nnz = " << matrixJA.size() << ") : matrix = " << storage )
        storage.prefetch();
        storageDense->assign( DenseStorage<ValueType>( numRows, numColumns, std::move( matrixDense ) ) );
        storage.wait();
        BOOST_CHECK_EQUAL( 0, storageDense->maxDiffNorm( storage ) );
        storage.purge();
        // Now set CSR data with an offset array
        utilskernel::HArrayUtils::scan1( matrixRowSizes, context );
        storage.setCSRData( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues );
        BOOST_CHECK_EQUAL( 0, storageDense->maxDiffNorm( storage ) );

        // ToDo: check for exceptions with illegal csr data, but only if SCAI_ASSERT_DEBUG is set
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( buildCSRSizesTest )
{
    using namespace hmemo;

    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    ContextPtr context = Context::getContextPtr();

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    IndexType numRows;
    IndexType numColumns;

    HArray<IndexType> matrixRowSizes;
    HArray<IndexType> matrixJA;
    HArray<ValueType> matrixValues;
    HArray<ValueType> matrixDense;

    getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );

    BOOST_REQUIRE_EQUAL( numRows, matrixRowSizes.size() );

    DenseStorage<ValueType> denseStorage( numRows, numColumns, std::move( matrixDense ) );

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_INFO( logger, "buildCSRSizes for " << storage )

        storage = denseStorage;   // convert DENSE -> current Format

        HArray<IndexType> rowSizes;

        storage.buildCSRSizes( rowSizes );

        BOOST_TEST( hostReadAccess( rowSizes ) == hostReadAccess( matrixRowSizes ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( buildCSRDataTest )
{
    using namespace hmemo;
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    ContextPtr context = Context::getContextPtr();
    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_INFO( logger, "buildCSRData for " << storage )
        storage.clear();
        IndexType numRows;
        IndexType numColumns;
        HArray<IndexType> matrixRowSizes;
        HArray<IndexType> matrixJA;
        HArray<ValueType> matrixValues;
        HArray<ValueType> matrixDense;
        getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
        IndexType numValues = matrixJA.size();
        // IA array not available yet, we have only sizes
        sparsekernel::CSRUtils::sizes2offsets( matrixRowSizes, matrixRowSizes, context );
        // Now we can use matrixRowSizes as IA array
        storage.setCSRData( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues );
        // make sure that we have all values stored
        BOOST_CHECK_EQUAL( numValues, storage.getNumValues() );
        SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                       << ", nnz = " << numValues << ") : storage = " << storage )

        HArray<IndexType> csrIA( context );
        HArray<IndexType> csrJA( context );
        HArray<ValueType> csrValues( context );

        storage.buildCSRData( csrIA, csrJA, csrValues );

        // check the IA array ( csrIA are offsets, matrixRowSizes was converted to offsets
 
        BOOST_TEST( hostReadAccess( matrixRowSizes ) == hostReadAccess( csrIA ), per_element() );
        BOOST_TEST( hostReadAccess( matrixJA ) == hostReadAccess( csrJA ), per_element() );
        BOOST_TEST( hostReadAccess( matrixValues ) == hostReadAccess( csrValues ), per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( fillCOOTest )
{
    using namespace hmemo;

    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    ContextPtr context = Context::getContextPtr();

    // input data (dense)

    const IndexType numRows    = 3;
    const IndexType numColumns = 4;

    HArray<ValueType> denseData( { 1, 0, 0, 4, 0, 2, 0, 0, 0, 0, 3, 0 } );

    // COO data to be filled

    HArray<IndexType> ia(     { 0, 1, 2, 2 } );
    HArray<IndexType> ja(     { 0, 2, 1, 1 } );
    HArray<ValueType> values( { 5, 6, 3, 3 } );

    // some COO arrays that should fail for filling, e.g. illegal size or indexes

    HArray<IndexType> badIA(     { 3, 1, 2, 2 } );
    HArray<IndexType> badJA(     { 0, 1, 5, 2 } );
    HArray<ValueType> badValues( { 0, 1, 5 } );        // illegal size

    // expected output data (dense)

    HArray<ValueType> denseCopy( { 5, 0, 0, 4, 0, 2, 6, 0, 0, 3, 3, 0 } );
    HArray<ValueType> denseAdd ( { 6, 0, 0, 4, 0, 2, 6, 0, 0, 6, 3, 0 } );

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    DenseStorage<ValueType> denseInput( numRows, numColumns, denseData );
    DenseStorage<ValueType> denseStorage;

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        if ( storage.getFormat() != Format::COO )
        {
            continue;  // REMOVE
        }

        storage.assign( denseInput );

        BOOST_CHECK_THROW(
        {
            storage.fillCOO( ia, ja, badValues, common::BinaryOp::COPY );
        }, Exception );

        BOOST_CHECK_THROW(
        {
            storage.fillCOO( badIA, ja, values, common::BinaryOp::COPY );
        }, Exception );

        BOOST_CHECK_THROW(
        {
            storage.fillCOO( ia, badJA, values, common::BinaryOp::COPY );
        }, Exception );

        storage.assign( denseInput );
        storage.fillCOO( ia, ja, values, common::BinaryOp::COPY );

        SCAI_LOG_DEBUG( logger, "storage filled (COPY) = " << storage )

        denseStorage.assign( storage );

        BOOST_TEST( hostReadAccess( denseStorage.getValues() ) == hostReadAccess( denseCopy ), 
                    boost::test_tools::per_element() );

        storage.assign( denseInput );
        storage.fillCOO( ia, ja, values, common::BinaryOp::ADD );

        SCAI_LOG_DEBUG( logger, "storage filled (ADD) = " << storage )

        denseStorage.assign( storage );

        BOOST_TEST( hostReadAccess( denseStorage.getValues() ) == hostReadAccess( denseAdd ), 
                    boost::test_tools::per_element() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
