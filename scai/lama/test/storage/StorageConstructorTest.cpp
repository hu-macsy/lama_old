/**
 * @file test/storage/StorageConstructorTest.cpp
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
 * @brief Contains constructor tests for storage classes
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( StorageConstructorTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.StorageConstructorTest" );

/* ------------------------------------------------------------------------- */

/** Define a list of all storagee types. */

typedef boost::mpl::list < 
            COOStorage<DefaultReal>,
            DIAStorage<DefaultReal>,
            CSRStorage<DefaultReal>,
            JDSStorage<DefaultReal>,
            ELLStorage<DefaultReal>,
            DenseStorage<DefaultReal>
        > StorageTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( defaultConstructorTest, StorageType, StorageTypes )
{
    StorageType storage;     // default constructor

    typedef typename StorageType::StorageValueType ValueType;
 
    // check zero sizes

    BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumRows() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumColumns() );

    // check correct format / type

    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, storage.getValueType() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sizeConstructorTest, StorageType, StorageTypes )
{
    hmemo::ContextPtr loc = hmemo::Context::getContextPtr();

    typedef typename StorageType::StorageValueType ValueType;

    const IndexType numRows = 13;
    const IndexType numCols = 17;

    // Attention: this constructor is DEPRECATED, use zero function

    StorageType storage( numRows, numCols, loc );

    BOOST_REQUIRE_EQUAL( numRows, storage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numCols, storage.getNumColumns() );

    ValueType expectedZero = 0;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            ValueType v = storage.getValue( i, j );
            BOOST_CHECK_EQUAL( expectedZero, v );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( zeroTest, StorageType, StorageTypes )
{
    hmemo::ContextPtr loc = hmemo::Context::getContextPtr();

    typedef typename StorageType::StorageValueType ValueType;

    const IndexType numRows = 13;
    const IndexType numCols = 17;

    StorageType storage = zero<StorageType>( numRows, numCols, loc );

    BOOST_REQUIRE_EQUAL( numRows, storage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numCols, storage.getNumColumns() );

    ValueType expectedZero = 0;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            auto v = storage.getValue( i, j );
            BOOST_CHECK_EQUAL( expectedZero, v );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, StorageType, StorageTypes )
{
    using namespace hmemo;

    ContextPtr context = Context::getContextPtr();

    typedef typename StorageType::StorageValueType ValueType;

    ValueType one = 1;

    auto storage1 = convert<StorageType>( DenseStorage<ValueType>( 2, 3, HArray<ValueType>( { 11, 12, 13, 21, 22, 23 } ) ) );
    auto storage2 = convert<StorageType>( DenseStorage<ValueType>( 3, 4, HArray<ValueType>( 12, one ) ) );

    storage1.swap( storage2 );

    BOOST_CHECK_EQUAL( 2, storage2.getNumRows() );
    BOOST_CHECK_EQUAL( 3, storage2.getNumColumns() );
    BOOST_CHECK_EQUAL( 3, storage1.getNumRows() );
    BOOST_CHECK_EQUAL( 4, storage1.getNumColumns() );

    BOOST_CHECK_EQUAL( 13, storage2.getValue( 0, 2 ) );
    BOOST_CHECK_EQUAL( 1, storage1.getValue( 0, 2 ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( newMatrixStorageTest, StorageType, StorageTypes )
{
    using namespace hmemo;

    const StorageType storage = zero<StorageType>( 3, 4 );

    std::unique_ptr<StorageType> storage1( storage.newMatrixStorage( 8, 7 ) );

    BOOST_CHECK_EQUAL( 8, storage1->getNumRows() );
    BOOST_CHECK_EQUAL( 7, storage1->getNumColumns() );

    const _MatrixStorage& s = storage;
    std::unique_ptr<_MatrixStorage> storage2( s.newMatrixStorage() );

    BOOST_CHECK_EQUAL( storage.getFormat(), storage2->getFormat() );
    BOOST_CHECK_EQUAL( storage.getNumRows(), storage2->getNumRows() );
    BOOST_CHECK_EQUAL( storage.getNumColumns(), storage2->getNumColumns() );

    std::unique_ptr<StorageType> storage3( storage.newMatrixStorage() );

    BOOST_CHECK_EQUAL( storage.getFormat(), storage3->getFormat() );
    BOOST_CHECK_EQUAL( storage.getNumRows(), storage3->getNumRows() );
    BOOST_CHECK_EQUAL( storage.getNumColumns(), storage3->getNumColumns() );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( TypeNameTest, StorageType, StorageTypes )
{   
    typedef typename StorageType::StorageValueType ValueType;

    std::string staticTypeName = StorageType::typeName();

    BOOST_CHECK( staticTypeName.length() > 0 );

    // typename, e.g. COOStorage<double>, should contain at least value type + format

    BOOST_CHECK( staticTypeName.find( common::TypeTraits<ValueType>::id() ) != std::string::npos );

    StorageType obj;

    std::string dynamicTypeName = obj.getTypeName();

    BOOST_CHECK_EQUAL( staticTypeName, dynamicTypeName );
    BOOST_CHECK( staticTypeName.find( format2Str( obj.getFormat() ) ) != std::string::npos );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
