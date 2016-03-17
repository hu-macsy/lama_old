/**
 * @file AllStorageTest.cpp
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
 * @brief Test cases applied to each storage class, i.e. test (virtual) methods of _MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/utilskernel/LArray.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/HArrayRef.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;
using utilskernel::LArray;

/* ------------------------------------------------------------------------------------------------------------------ */

/** This method checks for two arbitrary matrix storages if they are equal */

void checkEqual( const _MatrixStorage& storage1, const _MatrixStorage& storage2 )
{
    BOOST_REQUIRE_EQUAL( storage1.getNumRows(), storage2.getNumRows() );
    BOOST_REQUIRE_EQUAL( storage1.getNumColumns(), storage2.getNumColumns() );

    // take ScalarRepType that uses highest available precision

    LArray<ScalarRepType> row1;
    LArray<ScalarRepType> row2;

    for ( IndexType i = 0; i < storage1.getNumRows(); ++i )
    {
        storage1.getRow( row1, i );
        storage2.getRow( row2, i );

        hmemo::ReadAccess<ScalarRepType> rRow1( row1 );
        hmemo::ReadAccess<ScalarRepType> rRow2( row2 );

        for ( IndexType j = 0; j < storage1.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( rRow1[j], rRow2[j] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void initStorage( _MatrixStorage& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 5;

    static const double values[] =  { 6, 0, 7, 0, 0,
                                      0, 1, 0, 0, 0,
                                      0, 0, 9, 4, 0,
                                      2, 5, 0, 3, 8 };

    hmemo::HArrayRef<double> data( numRows * numColumns, values );
    DenseStorage<double> dense( data, numRows, numColumns );

    storage.assign( dense );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( AllStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.AllStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    Storages allMatrixStorages;    // is created by factory

    size_t nFormats = Format::UNDEFINED;
    size_t nTypes   = ARITHMETIC_HOST_TYPE_CNT;

    SCAI_LOG_INFO( logger, "Test all storages of factory to be empty, #storages = " << allMatrixStorages.size() )

    BOOST_CHECK_EQUAL( nTypes * nFormats, allMatrixStorages.size() );

    for ( size_t i = 0; i < allMatrixStorages.size(); ++i )
    {
        _MatrixStorage& storage = *allMatrixStorages[i];

        BOOST_CHECK_EQUAL( 0, storage.getNumRows() );
        BOOST_CHECK_EQUAL( 0, storage.getNumColumns() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for typeName" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];
    
        std::ostringstream s;

        s << storage;    // calls virtutal method writeAt for each storage class

        BOOST_CHECK( s.str().length() > 0 ); 
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( setIdentityTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    const IndexType N = 15;        // size of the matrix storage

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for setIdentity" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];

        storage.setIdentity(  N );

        SCAI_LOG_DEBUG( logger, "Identity matrix, N = " << N << " : " << storage )

        BOOST_REQUIRE_EQUAL( N, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumColumns() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumValues() );
        
        LArray<double> row;

        for ( IndexType i = 0; i < N; ++i )
        {
            storage.getRow( row, i );
            hmemo::ReadAccess<double> rRow( row );
        
            for ( IndexType j = 0; j < N; ++j )
            {
                if ( i == j )
                {
                    BOOST_CHECK_EQUAL( 1.0, rRow[j] );
                }
                else
                {
                    BOOST_CHECK_EQUAL( 0.0, rRow[j] );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( allocateTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    const IndexType numRows = 10;
    const IndexType numColumns = 15;

    ScalarRepType zero = 0;

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for setIdentity" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];

        storage.allocate( numRows, numColumns );

        SCAI_LOG_DEBUG( logger, "Zero matrix " << numRows << " x " << numColumns << " : " << storage )

        BOOST_REQUIRE_EQUAL( numRows, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( numColumns, storage.getNumColumns() );

        if ( storage.getFormat() == Format::DENSE )
        {
            // DenseStorage has all values and might be uninitialized

            continue;
        }
        
        BOOST_REQUIRE_EQUAL( 0, storage.getNumValues() );

        LArray<ScalarRepType> row;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            storage.getRow( row, i );
            hmemo::ReadAccess<ScalarRepType> rRow( row );
        
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( zero, rRow[j] );
            }
        }

        storage.clear();

        // verify that empty matrix has diagonal property

        BOOST_CHECK( storage.hasDiagonalProperty() );

        storage.allocate( 1, 1 );

        if ( storage.getFormat() == Format::DENSE )
        {
             // only dense matrix keeps its diagonal property
             BOOST_CHECK( storage.hasDiagonalProperty() );
         }
         else
         {
             BOOST_CHECK( !storage.hasDiagonalProperty() );
         }

         storage.purge();

         BOOST_CHECK( storage.hasDiagonalProperty() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( assignDenseTest )
{
    DenseStorage<double> denseStorage;

    initStorage( denseStorage );

    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages assign dense data" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];

        initStorage( storage );

        checkEqual( denseStorage, storage );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( conversionTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context

    Storages allMatrixStorages( context );    // is created by factory

    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for any conversion" )

    for ( size_t s1 = 0; s1 < allMatrixStorages.size(); ++s1 )
    {
        for ( size_t s2 = 0; s2 < allMatrixStorages.size(); ++s2 )
        {
            _MatrixStorage& storage1 = *allMatrixStorages[s1];
            _MatrixStorage& storage2 = *allMatrixStorages[s2];

            initStorage( storage1 );

            storage2 = storage1;   // converts both: type and format

            SCAI_LOG_DEBUG( logger, "Conversion: " << storage1 << " -> " << storage2 )

            checkEqual( storage1, storage2 );

            storage1.clear();

            // assignment of zero matrix, checks for consistent storage data

            storage2 = storage1;

            checkEqual( storage1, storage2 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END();
