/**
 * @file JDSStorageTest.cpp
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
 * @brief Test cases for JDSStorage( specific ones and all of MatrixStorageTest )
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/JDSStorage.hpp>
#include <lama/LAMAArrayUtils.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

namespace lama
{
namespace JDSStorageTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void checkTest( ContextPtr context )
{
    // This routine tests the check method of JDSStorage individually for this class
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
        const IndexType numValues = nJa;
        const IndexType numColumns = 10;
        const IndexType numDiagonals = nDlg;
        LAMAArrayRef<IndexType> jdsJA( nJa, valuesJa );
        LAMAArray<ValueType> jdsValues( nJa, 1.0 ); // only need to build JDS storage
        LAMAArrayRef<IndexType> jdsDLG( nDlg, valuesDlg );
        LAMAArrayRef<IndexType> jdsILG( nIlg, valuesIlg );
        LAMAArrayRef<IndexType> jdsPerm( nPerm, valuesPerm );
        JDSStorage<ValueType> jdsStorage;
        jdsStorage.setContext( context );
        // setJDSData will copy/convert values up to the needed context
        jdsStorage.setJDSData( numRows, numColumns, numValues, numDiagonals, jdsDLG, jdsILG, jdsPerm, jdsJA,
                               jdsValues );

        if ( icase == 0 )
        {
            jdsStorage.check( "test with correct values" );
        }
        else if ( icase == 1 )
        {
            //  -> invalid ja     { 0, 1, 2, 3, 15, 5, 6, 7, 8 }
            LAMAArray<IndexType>& jdsJA = const_cast<LAMAArray<IndexType>&>( jdsStorage.getJA() );
            LAMAArrayUtils::setVal( jdsJA, 5, 15 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal index in JA" ); }, Exception );
        }
        else if ( icase == 2 )
        {
            //  -> invalid ilg    { 3, 3, 4 }
            LAMAArray<IndexType>& jdsILG = const_cast<LAMAArray<IndexType>&>( jdsStorage.getIlg() );
            LAMAArrayUtils::setVal( jdsILG, 2, 4 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal ilg" ); }, Exception );
        }
        else if ( icase == 3 )
        {
            //  -> invalid dlg    { 3, 3, 4 }
            LAMAArray<IndexType>& jdsDLG = const_cast<LAMAArray<IndexType>&>( jdsStorage.getDlg() );
            LAMAArrayUtils::setVal( jdsDLG, 2, 4 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal dlg" ); }, Exception );
        }
        else if ( icase == 4 )
        {
            //  -> invalid perm   { 5, 0, 2 }
            LAMAArray<IndexType>& jdsPerm = const_cast<LAMAArray<IndexType>&>( jdsStorage.getPerm() );
            LAMAArrayUtils::setVal( jdsPerm, 0, 5 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal perm" ); }, Exception );
        }
        else if ( icase == 5 )
        {
            //  -> invalid perm   { 0, 0, 2 }
            LAMAArray<IndexType>& jdsPerm = const_cast<LAMAArray<IndexType>&>( jdsStorage.getPerm() );
            LAMAArrayUtils::setVal( jdsPerm, 0, 0 );
            BOOST_CHECK_THROW( { jdsStorage.check( "Expect illegal perm" ); }, Exception );
        }
    } // CASE_LOOP
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    JDSStorage<ValueType> jdsStorage;
    MatrixStorageTest<ValueType> storageTest( jdsStorage );
    storageTest.mMatrixStorage.setContext( loc );

    if ( base_test_case )
    {
        MATRIXSTORAGE_COMMONTESTCASES( storageTest );
    }
    else
    {
        storageTest.runTests();
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void constructorTest( )
{
    const IndexType numRows = 10;
    const IndexType numColumns = 15;
    JDSStorage<ValueType> jdsStorage( numRows, numColumns );
    BOOST_REQUIRE_EQUAL( numRows, jdsStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, jdsStorage.getNumColumns() );
    BOOST_REQUIRE_EQUAL( 0, jdsStorage.getNumValues() );

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

template<typename ValueType>
void constructorTest1(  ContextPtr context )
{
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
    LAMAArray<IndexType> jdsILG( numRows, ilg );
    LAMAArray<IndexType> jdsDLG( numDiagonals, dlg );
    LAMAArray<IndexType> jdsPerm( numRows, perm );
    LAMAArray<IndexType> jdsJA( numValues, ja );
    LAMAArray<ValueType> jdsValues( numValues, values );
    // Call the specific constructor for JDS storage
    JDSStorage<ValueType> jdsStorage( numRows, numColumns, numValues, numDiagonals,
                                      jdsDLG, jdsILG, jdsPerm, jdsJA, jdsValues );
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
    // copy constructor on all available locations
    JDSStorage<ValueType> jdsStorageCopy( jdsStorage, context );
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
void typeNameTest( )
{
    JDSStorage<ValueType> jdsStorage;
    std::string s = jdsStorage.typeName();
    BOOST_CHECK( s.length() > 0 );
}

} // namespace JDSStorageTest
} // namespace lama


/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( JDSStorageTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.JDSStorageTest" )

LAMA_AUTO_TEST_CASE_CT( checkTest, JDSStorageTest )
LAMA_AUTO_TEST_CASE_CT( commonTestCases, JDSStorageTest )
LAMA_AUTO_TEST_CASE_T( constructorTest, JDSStorageTest )
LAMA_AUTO_TEST_CASE_CT( constructorTest1, JDSStorageTest )
LAMA_AUTO_TEST_CASE_T( typeNameTest, JDSStorageTest )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
