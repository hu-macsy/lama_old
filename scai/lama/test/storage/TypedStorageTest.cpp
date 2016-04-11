/**
 * @file TypedStorageTest.cpp
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
 * @brief Test cases applied to each typed storage class, i.e. test (virtual) methods of MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/Context.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;
using utilskernel::LArray;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( TypedStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.TypedStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( factoryTest, ValueType, scai_arithmetic_test_types )
{
    TypedStorages<ValueType> allMatrixStorages;    // is created by factory

    size_t nFormats = Format::UNDEFINED;

    SCAI_LOG_INFO( logger, "factoryTest<" << common::TypeTraits<ValueType>::id() << "> : "
                   << allMatrixStorages.size() << " storages"                        )

    BOOST_CHECK_EQUAL( nFormats, allMatrixStorages.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( purgeTest, ValueType, scai_arithmetic_test_types )
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

BOOST_AUTO_TEST_CASE_TEMPLATE( normTest, ValueType, scai_arithmetic_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    SCAI_LOG_INFO( logger, "normTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        SCAI_LOG_DEBUG( logger, "normTest, storage = " << storage << " @ " << *storage.getContextPtr() )

        setDenseData( storage );

        ValueType maxNorm = storage.maxNorm();

        ValueType expected = 9.3f; // maximal absolute value

        SCAI_CHECK_CLOSE( maxNorm, expected, 1 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scaleTest, ValueType, scai_arithmetic_test_types )
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

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( diagonalTest, ValueType, scai_arithmetic_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        setDenseSquareData( storage );  // should be square

        SCAI_LOG_DEBUG( logger, "diagonalTest, storage = " << storage << " @ " << *storage.getContextPtr() )

        LArray<ValueType> diag( context );

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

        common::unique_ptr<MatrixStorage<ValueType> > inverse( storage.newMatrixStorage() );
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

        common::unique_ptr<MatrixStorage<ValueType> > inverse( storage.newMatrixStorage() );
        common::unique_ptr<MatrixStorage<ValueType> > result( storage.newMatrixStorage() );
        inverse->invert( storage );
        result->matrixTimesMatrix( ValueType( 1 ), storage, *inverse, ValueType( 0 ), *inverse );
        // BOOST_CHECK( result->maxDiffNorm( storage ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( inverseTestRandom1 )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        setDenseRandom( storage );

        common::unique_ptr<MatrixStorage<ValueType> > inverse( storage.copy() );
        common::unique_ptr<MatrixStorage<ValueType> > result( storage.newMatrixStorage() );
        inverse->invert( *inverse );
        result->matrixTimesMatrix( ValueType( 1 ), storage, *inverse, ValueType( 0 ), *inverse );
        // BOOST_CHECK( result->maxDiffNorm( storage ) < common::TypeTraits<ValueType>::small() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( matrixTimesVectorTest )
{
    // test of result = alpha * storage * x + beta * y
    // case 1: general
    // case 2: beta = 0, y is zero array
    // case 3: beta = 1, y == result, result += alpha * storage * x + beta

    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here

    const ValueType alpha = 1;
    const ValueType beta  = 2;
    const ValueType xVal  = 1.5;
    const ValueType yVal  = 0.5;

    DenseStorage<ValueType> denseStorage;

    setDenseData( denseStorage );

    const LArray<ValueType> x( denseStorage.getNumColumns(), xVal );
    const LArray<ValueType> y( denseStorage.getNumRows(), yVal );
    const LArray<ValueType> yDummy( 0, yVal );

    // for comparison: denseResult = alpha * DenseStorage * x + beta * y

    LArray<ValueType> denseResult1;
    LArray<ValueType> denseResult2;
    LArray<ValueType> denseResult3( y );

    denseStorage.matrixTimesVector( denseResult1, alpha, x, beta, y );
    denseStorage.matrixTimesVector( denseResult2, alpha, x, 0, yDummy );
    denseStorage.matrixTimesVector( denseResult3, alpha, x, 1, denseResult3 );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    SCAI_LOG_INFO( logger, "matrixTimesVectorTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        setDenseData( storage );

        SCAI_LOG_DEBUG( logger, "storage = " << storage )

        // result = alpha * storage * x + beta * y

        LArray<ValueType> result1( context );
        LArray<ValueType> result2( context );
        LArray<ValueType> result3( y, context );

        // BOOST_CHECK_THROW(
        // {
        //     storage.matrixTimesVector( result1, alpha, 1, beta, yDummy );
        // }, common::AssertException );

        storage.matrixTimesVector( result1, alpha, x, beta, y );
        storage.matrixTimesVector( result2, alpha, x, 0, yDummy );
        storage.matrixTimesVector( result3, alpha, x, 1, result3 );

        // should be the same as computed with dense storage

        BOOST_CHECK_EQUAL( denseResult1.maxDiffNorm( result1 ), 0 );
        BOOST_CHECK_EQUAL( denseResult2.maxDiffNorm( result2 ), 0 );
        BOOST_CHECK_EQUAL( denseResult3.maxDiffNorm( result3 ), 0 );
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

    setDenseData( denseStorage );

    const LArray<ValueType> x( denseStorage.getNumRows(), xVal );
    const LArray<ValueType> y( denseStorage.getNumColumns(), yVal );
    const LArray<ValueType> yDummy( 0, yVal );

    // for comparison: denseResult = alpha * x * DenseStorage + beta * y

    LArray<ValueType> denseResult1;
    LArray<ValueType> denseResult2;
    LArray<ValueType> denseResult3( y );

    denseStorage.vectorTimesMatrix( denseResult1, alpha, x, beta, y );
    denseStorage.vectorTimesMatrix( denseResult2, alpha, x, 0, yDummy );
    denseStorage.vectorTimesMatrix( denseResult3, alpha, x, 1, denseResult3 );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();

    SCAI_LOG_INFO( logger, "vectorTimesMatrixTest<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    TypedStorages<ValueType> allMatrixStorages( context );    // storage for each storage format

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];

        setDenseData( storage ); 

        SCAI_LOG_DEBUG( logger, "storage = " << storage )

        LArray<ValueType> result1( context );
        LArray<ValueType> result2( context );
        LArray<ValueType> result3( y, context );

        storage.vectorTimesMatrix( result1, alpha, x, beta, y );
        storage.vectorTimesMatrix( result2, alpha, x, 0, yDummy );
        storage.vectorTimesMatrix( result3, alpha, x, 1, result3 );

        // should be the same as computed with dense storage

        BOOST_CHECK_EQUAL( denseResult1.maxDiffNorm( result1 ), 0 );
        BOOST_CHECK_EQUAL( denseResult2.maxDiffNorm( result2 ), 0 );
        BOOST_CHECK_EQUAL( denseResult3.maxDiffNorm( result3 ), 0 );
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

        const LArray<ValueType> oldSolution( storage.getNumRows(), 1 );
        const LArray<ValueType> rhs( storage.getNumRows(), 2 );
        const LArray<ValueType> solution1;
        const LArray<ValueType> solution2;

        // clone the storage and set its diagonal to zero, but keep inverse of diagonal

        common::unique_ptr<MatrixStorage<ValueType> > storage1( storage.copy() );
        LArray<ValueType> diagonalInverse;
        storage1->getDiagonal( diagonalInverse );
        diagonalInverse.invert();
        storage1->setDiagonal( 0 );

        ValueType omegas[] = { 1.0, 0.8, 0.5 };
        const int NCASES = sizeof( omegas ) / sizeof( ValueType );

        for ( int k = 0; k < NCASES; ++k )
        {
            ValueType omega = omegas[k];

            LArray<ValueType> solution1( context );
            LArray<ValueType> solution2( context );

            storage.jacobiIterate( solution1, oldSolution, rhs, omega );

            const ValueType alpha = -1;
            const ValueType beta  = 1;

            //  solution2 = omega * ( rhs - B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution 

            storage1->matrixTimesVector( solution2, alpha, oldSolution, beta, rhs );
            solution2 *= diagonalInverse;
            utilskernel::HArrayUtils::arrayPlusArray( solution2, omega, solution2, 
                                                      ValueType( 1 ) - omega, oldSolution, solution2.getValidContext() );

            // solution1 and solution2 must be the same
    
            BOOST_CHECK( solution1.maxDiffNorm( solution2 ) < common::TypeTraits<ValueType>::small() );
        }
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

        SCAI_LOG_INFO( logger, "setCSRData for " << storage )

        storage.clear();
        IndexType numRows;
        IndexType numColumns;
        LArray<IndexType> matrixRowSizes( context );
        LArray<IndexType> matrixJA( context );
        LArray<ValueType> matrixValues( context );
        LArray<ValueType> matrixDense( context );
        getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
        IndexType numValues = matrixJA.size();
        _MatrixStorage::sizes2offsets( matrixRowSizes );

        // Now we can use matrixRowSizes as IA array

        storage.setCSRData( numRows, numColumns, numValues, matrixRowSizes, matrixJA, matrixValues );

        SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                       << ", nnz = " << numValues << ") : matrix = " << storage )
        storage.prefetch();
        hmemo::ReadAccess<ValueType> results( matrixDense );

        for ( IndexType i = 0; i < numRows; i++ )
        {
            for ( IndexType j = 0; j < numColumns; j++ )
            {
                BOOST_CHECK_EQUAL( storage.getValue( i, j ) , results[ i * numColumns + j ] );
            }
        }

        storage.wait();

        // now we make some checks with incorrect CSR data

        /* not yet

        BOOST_CHECK_THROW(
        {
            storage.setCSRData( numRows, numColumns, numValues-1, matrixRowSizes, matrixJA, matrixValues );
        },
        common::Exception );

        */

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
        LArray<IndexType> matrixRowSizes;
        LArray<IndexType> matrixJA;
        LArray<ValueType> matrixValues;
        LArray<ValueType> matrixDense;
        getMatrix_7_4 ( numRows, numColumns, matrixRowSizes, matrixJA, matrixValues, matrixDense );
        IndexType numValues = matrixJA.size();
        // IA array not available yet, we have only sizes
        _MatrixStorage::sizes2offsets( matrixRowSizes );
        // Now we can use matrixRowSizes as IA array
        storage.setCSRData( numRows, numColumns, numValues, matrixRowSizes, matrixJA, matrixValues );
        // make sure that storage has not diagonal property, otherwise it will build wrong CSR data
        BOOST_REQUIRE_EQUAL( storage.hasDiagonalProperty(), false );
        // make sure that we have all values stored
        BOOST_CHECK_EQUAL( numValues, storage.getNumValues() );

        SCAI_LOG_INFO( logger, "set CSR data (" << numRows << " x " << numColumns
                       << ", nnz = " << numValues << ") : storage = " << storage )

        LArray<IndexType> csrIA( context );
        LArray<IndexType> csrJA( context );
        LArray<ValueType> csrValues( context );

        storage.buildCSRData( csrIA, csrJA, csrValues );

        BOOST_REQUIRE_EQUAL( numRows + 1, csrIA.size() );
        BOOST_REQUIRE_EQUAL( numValues, csrJA.size() );
        BOOST_REQUIRE_EQUAL( numValues, csrValues.size() );

        // check the IA array ( csrIA are offsets, matrixRowSizes was converted to offsets

        BOOST_CHECK_EQUAL( 0, matrixRowSizes.maxDiffNorm( csrIA ) );
        BOOST_CHECK_EQUAL( 0, matrixJA.maxDiffNorm( csrJA ) );
        BOOST_CHECK_EQUAL( 0, matrixValues.maxDiffNorm( csrValues ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
