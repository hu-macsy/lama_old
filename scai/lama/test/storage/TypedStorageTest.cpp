/**
 * @file TypedStorageTest.cpp
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
 * @brief Test cases applied to each typed storage class, i.e. test (virtual) methods of MatrixStorage
 * @author Thomas Brandes
 * @date 31.08.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/hmemo/ReadAccess.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( TypedStorageTest );

SCAI_LOG_DEF_LOGGER( logger, "Test.TypedStorageTest" )

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_CASE_TEMPLATE( factoryTest, ValueType, scai_arithmetic_test_types )
{
    TypedStorages<ValueType> allMatrixStorages;    // is created by factory

    size_t nFormats = Format::UNDEFINED;

    SCAI_LOG_INFO( logger, "factoryTest<" << common::TypeTraits<ValueType>::id() << "> : " 
                            << allMatrixStorages.size() << " storages"                        )

    BOOST_CHECK_EQUAL( nFormats, allMatrixStorages.size() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
