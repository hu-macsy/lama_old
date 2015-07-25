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
 * @brief Contains the implementation of the class DenseStorageTest
 * @author Alexander Büchel
 * @date 12.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/storage/DenseStorage.hpp>

#include <test/MatrixStorageTest.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

extern bool base_test_case;
extern std::string testcase;

namespace lama
{
namespace DenseStorageTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    DenseStorage<ValueType> denseStorage;
    MatrixStorageTest<ValueType> storageTest( denseStorage );
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
void typeNameTest()
{
    DenseStorage<ValueType> denseStorage;
    std::string s = denseStorage.typeName();
    BOOST_CHECK( s.length() > 0 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setZeroTest()
{
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

} // namespace DenseStorageTest
} // namespace lama

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( DenseStorageTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.DenseStorageTest" )
LAMA_AUTO_TEST_CASE_CT( commonTestCases, DenseStorageTest )
LAMA_AUTO_TEST_CASE_T( setZeroTest, DenseStorageTest )
LAMA_AUTO_TEST_CASE_T( typeNameTest, DenseStorageTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
