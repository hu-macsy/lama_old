/**
 * @file AllStorageTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Test cases applied to each storage class, i.e. test (virtual) methods of _MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/HArrayRef.hpp>

#include <scai/logging.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/macros/count.hpp>

using namespace scai;
using namespace lama;
using hmemo::HArray;

using boost::test_tools::per_element;

/* ------------------------------------------------------------------------------------------------------------------ */

/** This method checks for two arbitrary matrix storages if they are equal */

void checkEqual( const _MatrixStorage& storage1, const _MatrixStorage& storage2 )
{
    using namespace hmemo;

    auto dense1 = convert<DenseStorage<ScalarRepType>>( storage1 );
    auto dense2 = convert<DenseStorage<ScalarRepType>>( storage2 );

    BOOST_TEST( hostReadAccess( dense1.getValues() ) == hostReadAccess( dense2.getValues() ), per_element() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

void initStorage( _MatrixStorage& storage )
{
    typedef SCAI_TEST_TYPE ValueType;
    const IndexType numRows = 4;
    const IndexType numColumns = 5;
    static const ValueType values[] =  { 6, 0, 7, 0, 0,
                                         0, 1, 0, 0, 0,
                                         0, 0, 9, 4, 0,
                                         2, 5, 0, 3, 8
                                       };
    hmemo::HArrayRef<ValueType> data( numRows * numColumns, values );
    DenseStorage<ValueType> dense( numRows, numColumns, data );
    storage.assign( dense );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( AllStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.AllStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    Storages allMatrixStorages;    // is created by factory
    size_t nFormats = static_cast<size_t>( Format::UNDEFINED );
    nFormats -= 1; // not for STENCIL
    size_t nTypes   = SCAI_COMMON_COUNT_NARG( SCAI_NUMERIC_TYPES_HOST );
    SCAI_LOG_INFO( logger, "#formats = " << nFormats << ", #types = " << nTypes )
    SCAI_LOG_INFO( logger, "Test all storages of factory to be empty, #storages = " << allMatrixStorages.size() )
    BOOST_CHECK_EQUAL( nTypes * nFormats, allMatrixStorages.size() );

    for ( size_t i = 0; i < allMatrixStorages.size(); ++i )
    {
        _MatrixStorage& storage = *allMatrixStorages[i];
        BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumRows() );
        BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumColumns() );
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
        std::ostringstream os;
        os << storage;    // calls virtutal method writeAt for each storage class
        BOOST_CHECK( os.str().length() > 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( setIdentityTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    const IndexType N = 15;        // size of the matrix storage
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for setIdentity" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        storage.setIdentity(  N );
        SCAI_LOG_DEBUG( logger, "Identity matrix, N = " << N << " : " << storage )
        BOOST_REQUIRE_EQUAL( N, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumColumns() );
        BOOST_REQUIRE_EQUAL( N, storage.getNumValues() );
        HArray<ValueType> row;

        for ( IndexType i = 0; i < N; ++i )
        {
            storage.getRow( row, i );
            hmemo::ReadAccess<ValueType> rRow( row );

            for ( IndexType j = 0; j < N; ++j )
            {
                if ( i == j )
                {
                    BOOST_CHECK_EQUAL( ValueType( 1 ), rRow[j] );
                }
                else
                {
                    BOOST_CHECK_EQUAL( ValueType( 0 ), rRow[j] );
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( allocateTest, ValueType, scai_numeric_test_types )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    ValueType zero = 0;
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for setIdentity" )

    HArray<ValueType> zeroColumn( numColumns, zero );

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        storage.allocate( numRows, numColumns );
        SCAI_LOG_DEBUG( logger, "Zero matrix " << numRows << " x " << numColumns << " : " << storage )
        BOOST_REQUIRE_EQUAL( numRows, storage.getNumRows() );
        BOOST_REQUIRE_EQUAL( numColumns, storage.getNumColumns() );
        BOOST_REQUIRE_EQUAL( IndexType( 0 ), storage.getNumValues() );
        HArray<ValueType> row;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            storage.getRow( row, i );

            BOOST_TEST( hostReadAccess( row ) == hostReadAccess( zeroColumn ), per_element() );
        }

        storage.clear();
        storage.purge();
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( assignDenseTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    DenseStorage<ValueType> denseStorage;
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

BOOST_AUTO_TEST_CASE_TEMPLATE( getColTest, ValueType, scai_numeric_test_types )
{
    DenseStorage<ValueType> denseStorage;
    initStorage( denseStorage );
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    TypedStorages<ValueType> allMatrixStorages( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages assign dense data" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        initStorage( storage );

        HArray<ValueType> col1;
        HArray<ValueType> col2;

        for ( IndexType j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            storage.getColumn( col1, j );
            denseStorage.getColumn( col2, j );

            BOOST_TEST( hostReadAccess( col1 ) == hostReadAccess( col2 ), per_element() );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( printTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Storages allMatrixStorages( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages assign dense data" )

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        _MatrixStorage& storage = *allMatrixStorages[s];
        initStorage( storage );
        std::ostringstream os;
        storage.print( os );
        BOOST_CHECK( os.str().length()  > 0 );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( conversionTest )
{
    hmemo::ContextPtr context = hmemo::Context::getContextPtr();  // test context
    Storages allMatrixStorages( context );    // is created by factory
    SCAI_LOG_INFO( logger, "Test " << allMatrixStorages.size() << "  storages for any conversion" )

    const size_t n = allMatrixStorages.size();

    for ( size_t s1 = 0; s1 < n; ++s1 )
    {
        for ( size_t s2 = 0; s2 < n; ++s2 )
        {
            _MatrixStorage& storage1 = *allMatrixStorages[s1];
            _MatrixStorage& storage2 = *allMatrixStorages[s2];

            initStorage( storage1 );

            SCAI_LOG_INFO( logger, "conversionTest " << s1 << " x " << s2 << " of " << n << " x " << n 
                                   << ", storage1 = " << storage1 << ", -> storage2  = " << storage2 )


            storage2 = storage1;   // converts both: type and format
            checkEqual( storage1, storage2 );

            storage1.clear();

            // assignment of zero matrix, checks for consistent storage data

            storage2 = storage1;
            checkEqual( storage1, storage2 );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
