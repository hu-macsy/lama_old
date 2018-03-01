/**
 * @file DenseStorageTest.cpp
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
 * @brief Test cases for DenseStorage( only specific ones )
 * @author Thomas Brandes
 * @date 12.03.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( DenseStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseStorageTest" )

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( setZeroTest, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "setZeroTest for DenseStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] =
    {
        6.0, 0.0, 0.0, 4.0,
        7.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 9.0, 4.0,
        2.0, 5.0, 0.0, 3.0
    };
    DenseStorage<ValueType> denseStorage;
    denseStorage.setContextPtr( context );
    denseStorage.setRawDenseData( numRows, numColumns, values );
    denseStorage.setZero();

    for ( IndexType i = 0; i < denseStorage.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( denseStorage.getValue( i, j ), 0.0 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "constructorTest for DenseStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    const IndexType numRows = 4;
    const IndexType numColumns = 2;
    // DenseStorage<ValueType> denseStorage( numRows, numColumns );
    DenseStorage<ValueType> denseStorage;
    denseStorage.setContextPtr( context );
    denseStorage.allocate( numRows, numColumns );
    denseStorage.setZero();
    BOOST_REQUIRE_EQUAL( numRows, denseStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, denseStorage.getNumColumns() );

    for ( IndexType i = 0; i < denseStorage.getNumRows(); ++i )
    {
        for ( IndexType j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( denseStorage.getValue( i, j ), ValueType( 0 ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
static inline const ValueType* getPointer( const HArray<ValueType>& a, ContextPtr ctx )
{
    ReadAccess<ValueType> rA( a, ctx );
    return rA.get();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( splitUpTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    ContextPtr context = Context::getContextPtr();

    const IndexType numRows    = 4;
    const IndexType numColumns = 3;

    HArray<ValueType> denseValues( { 5, 0, 0, 0, 5, 0, 0, 3, 2, 0, 1, 1 } );

    const ValueType* ptrValues = getPointer( denseValues, context );

    DenseStorage<ValueType> denseStorage( numRows, numColumns, std::move( denseValues ) );

    BOOST_CHECK_EQUAL( ptrValues, getPointer( denseStorage.getData(), context ) );

    IndexType outNumRows;
    IndexType outNumColumns;
    HArray<ValueType> outValues;

    denseStorage.splitUp( outNumRows, outNumColumns, outValues );

    BOOST_CHECK_EQUAL( outNumRows, numRows );
    BOOST_CHECK_EQUAL( outNumColumns, numColumns );

    // we have still the same dense data as before moving it 

    BOOST_CHECK_EQUAL( ptrValues, getPointer( outValues, context ) );

    // verify that the remaining storage is nothing else than a dense storge

    DenseStorage<ValueType> zeroDenseStorage;
    BOOST_CHECK_EQUAL( zeroDenseStorage.maxDiffNorm( denseStorage ), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for DenseStorage<" << common::TypeTraits<ValueType>::id() << ">" )
    storageTypeNameTest<DenseStorage<ValueType> >( "Dense" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( DenseCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<DenseStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
