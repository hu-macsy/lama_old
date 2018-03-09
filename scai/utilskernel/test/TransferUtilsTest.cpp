/**
 * @file HArrayUtilsTest.cpp
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
 * @brief Tests for the class TransferUtils
 * @author Thomas Brandes
 * @date 22.01.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/TransferUtils.hpp>
#include <scai/utilskernel.hpp>

#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/utilskernel/test/HArrays.hpp>

#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>

#include <typeinfo>
#include <memory>

using std::unique_ptr;

using namespace scai;
using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TransferUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.TransferUtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    auto ctx = Context::getContextPtr();

    HArray<ValueType> sourceArray( { 1, 2, 3, 4, 5, 6 }, ctx );
    HArray<ValueType> targetArray( { 9, 9, 9, 9, 9 } );

    HArray<IndexType> sourceIndexes( { 2, 3, 1, 5 } );
    HArray<IndexType> targetIndexes( { 4, 2, 1, 0 } );

    SCAI_ASSERT_EQ_ERROR( sourceIndexes.size(), targetIndexes.size(), "#indexes must match" )

    HArray<ValueType> expectedTargetArray( { 6, 2, 4, 9, 3 } );

    TransferUtils::copy( targetArray, targetIndexes, sourceArray, sourceIndexes );

    BOOST_TEST( hostReadAccess( expectedTargetArray ) == hostReadAccess( targetArray ),
                 boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( gatherVTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    IndexType nRows = 5;

    const IndexType rowOffsets[] = { 0,    2,       5, 6, 7,    9 };
    const ValueType vals[]       = { 2, 1, 3, 1, 5, 6, 1, 2, 3 };

    const IndexType sIndexes[]   = { 3, 1 };

    const IndexType nTarget      = 4;
    const ValueType eVals[]      = { 1, 3, 1, 5 };

    HArray<ValueType> sourceArray( rowOffsets[nRows], vals );
    HArray<IndexType> sourceOffsets( nRows + 1, rowOffsets );
    HArray<IndexType> sourceIndexes( 2, sIndexes );

    HArray<ValueType> targetArray( 4 );

    TransferUtils::gatherV( targetArray, sourceArray, sourceOffsets, sourceIndexes );

    ReadAccess<ValueType> targetRead( targetArray );

    for ( IndexType i = 0; i < nTarget; ++i ) 
    {
        BOOST_CHECK_EQUAL( targetRead[i], eVals[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

