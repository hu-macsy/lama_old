/**
 * @file DenseStorageTest.cpp
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
 * @brief Test cases for DenseStorage( only specific ones )
 * @author Thomas Brandes
 * @date 12.03.2012
 * @since 1.0.0
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

BOOST_AUTO_TEST_CASE_TEMPLATE( setZeroTest, ValueType, scai_arithmetic_test_types )
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
 
    ValueType eps = static_cast<ValueType>( 1E-5 );

    DenseStorage<ValueType> denseStorage;
    denseStorage.setContextPtr( context );

    denseStorage.setRawDenseData( numRows, numColumns, values, eps );
    denseStorage.setZero();

    for ( int i = 0; i < denseStorage.getNumRows(); ++i )
    {
        for ( int j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( denseStorage.getValue( i, j ), 0.0 );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_arithmetic_test_types )
{
    ContextPtr context = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "constructorTest for DenseStorage<" << common::TypeTraits<ValueType>::id() << "> @ " << *context )

    const IndexType numRows = 4;
    const IndexType numColumns = 2;

    DenseStorage<ValueType> denseStorage( numRows, numColumns );
    denseStorage.setContextPtr( context );
    denseStorage.setZero();

    BOOST_REQUIRE_EQUAL( numRows, denseStorage.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, denseStorage.getNumColumns() );

    for ( int i = 0; i < denseStorage.getNumRows(); ++i )
    {
        for ( int j = 0; j < denseStorage.getNumColumns(); ++j )
        {
            BOOST_CHECK_EQUAL( denseStorage.getValue( i, j ), ValueType( 0 ) );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_arithmetic_test_types )
{
    SCAI_LOG_INFO( logger, "swapTest for DenseStorage<" << common::TypeTraits<ValueType>::id() << ">" )

    // use template storage test 

    storageSwapTest<DenseStorage<ValueType> >();
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( typenameTest, ValueType, scai_arithmetic_test_types )
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
