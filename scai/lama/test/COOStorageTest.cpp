/**
 * @file COOStorageTest.cpp
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
 * @brief Contains the implementation of the class COOStorageTest
 * @author Alexander Büchel
 * @date 12.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/storage/COOStorage.hpp>

#include <scai/lama/test/MatrixStorageTest.hpp>
#include <scai/lama/test/TestMacros.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>

using namespace scai::utilskernel;
using namespace scai::lama;
using namespace scai::hmemo;
using scai::common::Exception;

extern bool base_test_case;
extern std::string testcase;

namespace scai
{
namespace lama
{
namespace COOStorageTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    COOStorage<ValueType> cooStorage;
    MatrixStorageTest<ValueType> storageTest( cooStorage );
    storageTest.mMatrixStorage.setContextPtr( loc );

    if ( base_test_case )
    {
        MATRIXSTORAGE_COMMONTESTCASES( storageTest );
    }
    else
    {
        storageTest.runTests();
    }
}

} /* end namespace COOStorageTest */
} /* end namespace lama */
} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( COOStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.COOStorageTest" )
LAMA_AUTO_TEST_CASE_CT( commonTestCases, COOStorageTest, scai::lama )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
