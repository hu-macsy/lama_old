/**
 * @file JDSStorageTest.cpp
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
 * @brief Test cases for JDSStorage( only specific ones )
 * @author Thomas Brandes
 * @date 31.08.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/JDSStorage.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/logging.hpp>

#include <scai/lama/test/storage/StorageTemplateTests.hpp>

using namespace scai;
using namespace scai::lama;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using scai::common::Exception;
using scai::common::BinaryOp;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( JDSStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.JDSStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( checkTest, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    // This routine tests the check method of JDSStorage individually for this class
    JDSStorage<ValueType> jdsStorage;
    jdsStorage.setContextPtr( context );
    SCAI_LOG_INFO( logger, "checkTest for JDSStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    for ( int icase = 0; icase < 6; ++icase )
    {
        // build up a correct JDS storage
        IndexType valuesJa[] =
        { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
        const IndexType nJa = sizeof( valuesJa ) / sizeof( IndexType );
        IndexType valuesDlg[] =
        { 3, 3, 3 };
        const IndexType nDlg = sizeof( valuesDlg ) / sizeof( IndexType );
        IndexType valuesIlg[] =
        { 3, 3, 3 };
        const IndexType nIlg = sizeof( valuesIlg ) / sizeof( IndexType );
        IndexType valuesPerm[] =
        { 2, 0, 1 };
        const IndexType nPerm = sizeof( valuesPerm ) / sizeof( IndexType );
        const IndexType numRows = nIlg;
        const IndexType numColumns = 10;
        HArrayRef<IndexType> jdsJA( nJa, valuesJa );
        HArray<ValueType> jdsValues( nJa, 1.0 ); // only need to build JDS storage
        HArrayRef<IndexType> jdsDLG( nDlg, valuesDlg );
        HArrayRef<IndexType> jdsILG( nIlg, valuesIlg );
        HArrayRef<IndexType> jdsPerm( nPerm, valuesPerm );
        // setJDSData will copy/convert values up to the needed context

        jdsStorage.setJDSData( numRows, numColumns, jdsDLG, jdsILG, jdsPerm, jdsJA, jdsValues );

        if ( icase == 0 )
        {
            jdsStorage.check( "test with correct values" );
        }
        else if ( icase == 1 )
        {
            //  -> invalid ja     { 0, 1, 2, 3, 15, 5, 6, 7, 8 }
            HArray<IndexType>& jdsJA = const_cast<HArray<IndexType>&>( jdsStorage.getJA() );
            IndexType idx = 5;
            IndexType val = 15;
            HArrayUtils::setVal( jdsJA, idx, val );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ilg    { 3, 3, 4 }
            HArray<IndexType>& jdsILG = const_cast<HArray<IndexType>&>( jdsStorage.getIlg() );
            IndexType idx = 2;
            IndexType val = 4;
            HArrayUtils::setVal( jdsILG, idx, val, BinaryOp::COPY );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal ilg" ); }, Exception );
        }
        else if ( icase == 3 )
        {
            //  -> invalid dlg    { 3, 3, 4 }
            HArray<IndexType>& jdsDLG = const_cast<HArray<IndexType>&>( jdsStorage.getDlg() );
            IndexType idx = 2;
            IndexType val = 4;
            HArrayUtils::setVal( jdsDLG, idx, val, BinaryOp::COPY );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal dlg" ); }, Exception );
        }
        else if ( icase == 4 )
        {
            //  -> invalid perm   { 5, 0, 2 }
            HArray<IndexType>& jdsPerm = const_cast<HArray<IndexType>&>( jdsStorage.getPerm() );
            HArrayUtils::setVal<IndexType>( jdsPerm, 0, 5 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal perm" ); }, Exception );
        }
        else if ( icase == 5 )
        {
            //  -> invalid perm   { 0, 0, 2 }
            HArray<IndexType>& jdsPerm = const_cast<HArray<IndexType>&>( jdsStorage.getPerm() );
            IndexType idx = 0;
            IndexType val = 0;
            HArrayUtils::setVal( jdsPerm, idx, val, BinaryOp::COPY );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal perm" ); }, Exception );
        }
    } // CASE_LOOP
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "constructorTest for JDSStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    JDSStorage<ValueType> jdsStorage;
    jdsStorage.setContextPtr( context );
    jdsStorage.allocate( numRows, numColumns );
    BOOST_REQUIRE_EQUAL( numRows, jdsStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( IndexType( 0 ), jdsStorage.getNumValues() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            float v = static_cast<float> ( jdsStorage.getValue( i, j ) );
            BOOST_CHECK_SMALL( v , 1.0e-5f );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructor1Test, ValueType, scai_numeric_test_types )
{
    ContextPtr context = Context::getContextPtr();
    SCAI_LOG_INFO( logger, "constructor1Test for JDSStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )
    const IndexType numRows = 3;
    const IndexType numColumns = 3;
// Note: ja, values are stored column-major order
    const IndexType ja[] =
    {   2, 0, 1, 2,};
    const ValueType values[] =
    {   0.3f, 0.5f, 0.5f, 0.2f};
    const IndexType perm[] =
    {   2, 0, 1};
    const IndexType ilg[] =
    {   2, 1, 1};
    const IndexType dlg[] =
    {   3, 1};
    const IndexType numDiagonals = sizeof( dlg ) / sizeof( IndexType );
    const IndexType numValues = sizeof( ja ) / sizeof( IndexType );
    const IndexType sizeValues = sizeof( values ) / sizeof( ValueType );
    const IndexType sizePerm = sizeof( perm ) / sizeof( IndexType );
    const IndexType sizeILG = sizeof( ilg ) / sizeof( IndexType );
    BOOST_CHECK_EQUAL( numValues, sizeValues );
    BOOST_CHECK_EQUAL( numRows, sizePerm );
    BOOST_CHECK_EQUAL( numRows, sizeILG );
    HArray<IndexType> jdsILG( numRows, ilg );
    HArray<IndexType> jdsDLG( numDiagonals, dlg );
    HArray<IndexType> jdsPerm( numRows, perm );
    HArray<IndexType> jdsJA( numValues, ja );
    HArray<ValueType> jdsValues( numValues, values );

    // Call the specific constructor for JDS storage

    JDSStorage<ValueType> jdsStorage( numRows, numColumns, 
                                      jdsDLG, jdsILG, jdsPerm, jdsJA, jdsValues, context );

    // Test all the getter routines, includes the specific ones for JDSStorage
    BOOST_REQUIRE_EQUAL( numRows, jdsStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numDiagonals, jdsStorage.getNumDiagonals() );
    {
        ReadAccess<IndexType> jdsILG( jdsStorage.getIlg() );
        ReadAccess<IndexType> jdsPerm( jdsStorage.getPerm() );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( ilg[i], jdsILG[i] );
            BOOST_CHECK_EQUAL( perm[i], jdsPerm[i] );
        }

        ReadAccess<IndexType> jdsDLG( jdsStorage.getDlg() );

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            BOOST_CHECK_EQUAL( dlg[i], jdsDLG[i] );
        }

        ReadAccess<IndexType> jdsJA( jdsStorage.getJA() );
        ReadAccess<ValueType> jdsValues( jdsStorage.getValues() );

        // JDS keeps values in same order

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], jdsJA[i] );
            BOOST_CHECK_EQUAL( values[i], jdsValues[i] );
        }
    }

    // copy constructor 

    JDSStorage<ValueType> jdsStorageCopy( jdsStorage );
    BOOST_REQUIRE_EQUAL( numRows, jdsStorageCopy.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorageCopy.getNumColumns() );
    BOOST_REQUIRE_EQUAL( numDiagonals, jdsStorageCopy.getNumDiagonals() );
    {
        ReadAccess<IndexType> jdsILG( jdsStorageCopy.getIlg() );
        ReadAccess<IndexType> jdsPerm( jdsStorageCopy.getPerm() );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            BOOST_CHECK_EQUAL( ilg[i], jdsILG[i] );
            BOOST_CHECK_EQUAL( perm[i], jdsPerm[i] );
        }

        ReadAccess<IndexType> jdsDLG( jdsStorageCopy.getDlg() );

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            BOOST_CHECK_EQUAL( dlg[i], jdsDLG[i] );
        }

        ReadAccess<IndexType> jdsJA( jdsStorageCopy.getJA() );
        ReadAccess<ValueType> jdsValues( jdsStorageCopy.getValues() );

        // JDS keeps values in same order

        for ( IndexType i = 0; i < numValues; ++i )
        {
            BOOST_CHECK_EQUAL( ja[i], jdsJA[i] );
            BOOST_CHECK_EQUAL( values[i], jdsValues[i] );
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

    const IndexType numRows = 3;
    const IndexType numColumns = 3;
    const IndexType numDiagonals = 2;
    const IndexType numValues = 4;

    HArray<IndexType> jdsIlg( { 2, 1, 1 } );
    HArray<IndexType> jdsDlg( { 3, 1 } );
    HArray<IndexType> jdsPerm( { 2, 0, 1 } );
    HArray<IndexType> jdsJA( { 2, 0, 1, 2 } );
    HArray<ValueType> jdsValues( { 3, 5, 5, 2 } );

    const IndexType* ptrJA = getPointer( jdsJA, context );
    const ValueType* ptrValues = getPointer( jdsValues, context );

    JDSStorage<ValueType> jdsStorage( numRows, numColumns, 
                                      std::move( jdsDlg ), std::move( jdsIlg ), std::move( jdsPerm ),
                                      std::move( jdsJA ), std::move( jdsValues ) );

    IndexType outNumRows;
    IndexType outNumColumns;

    HArray<IndexType> outIlg;
    HArray<IndexType> outDlg;
    HArray<IndexType> outPerm;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;

    jdsStorage.splitUp( outNumRows, outNumColumns, outDlg, outIlg, outPerm, outJA, outValues );

    BOOST_CHECK_EQUAL( outDlg.size(), numDiagonals );
    BOOST_CHECK_EQUAL( outIlg.size(), numRows );
    BOOST_CHECK_EQUAL( outPerm.size(), numRows );
    BOOST_CHECK_EQUAL( outJA .size(), numValues );
    BOOST_CHECK_EQUAL( outValues.size(), numValues );

    BOOST_CHECK_EQUAL( outNumRows, numRows );
    BOOST_CHECK_EQUAL( outNumColumns, numColumns );

    BOOST_CHECK_EQUAL( ptrJA, getPointer( outJA, context ) );
    BOOST_CHECK_EQUAL( ptrValues, getPointer( outValues, context ) );

    BOOST_CHECK_EQUAL( jdsStorage.getIlg().size(), 0 );
    BOOST_CHECK_EQUAL( jdsStorage.getDlg().size(), 0 );
    BOOST_CHECK_EQUAL( jdsStorage.getJA().size(), 0 );
    BOOST_CHECK_EQUAL( jdsStorage.getValues().size(), 0 );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_numeric_test_types )
{
    SCAI_LOG_INFO( logger, "typeNameTest for JDSStorage<" << common::TypeTraits<ValueType>::id() << ">" )
    storageTypeNameTest<JDSStorage<ValueType> >( "JDS" );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE( JDSCopyTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // test for one value type is sufficient here
    copyStorageTest<JDSStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
