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

#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/lama/test/MatrixStorageTest.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using scai::common::Exception;

extern bool base_test_case;
extern std::string testcase;

namespace scai
{
namespace lama
{
namespace JDSStorageTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    JDSStorage<ValueType> jdsStorage;
    MatrixStorageTest<ValueType> storageTest( jdsStorage );
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

} /* end namespace JDSStorageTest */

} /* end namespace lama */

} /* end namespace scai */


/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( JDSStorageTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.JDSStorageTest" )

LAMA_AUTO_TEST_CASE_CT( commonTestCases, JDSStorageTest, scai::lama )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
