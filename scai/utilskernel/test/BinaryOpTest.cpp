/**
 * @file BinaryOpTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Test enum for BinaryOp
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <sstream>

using namespace scai;
using namespace utilskernel;

BOOST_AUTO_TEST_CASE( BinaryOpTest )
{
    for ( int type = binary::COPY; type < binary::MAX_BINARY_OP; ++type )
    {
        std::ostringstream s;
        s << binary::BinaryOp( type );
        BOOST_CHECK( s.str().length() > 0 );

        if ( type == binary::COPY )
        {
            BOOST_CHECK_EQUAL( s.str(), "COPY" );
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( ApplyBinOpTest, T, scai_array_test_types )
{
    // we do not check here the results but make sure that all ops are supported

    for ( int op = binary::COPY; op < binary::MAX_BINARY_OP; ++op )
    {
        T x = 5;
        T y = 2;
        T z = 0;

        if ( common::TypeTraits<T>::stype == common::TypeTraits<IndexType>::stype )
        {
            if ( op == binary::POW || op == binary::COPY_SIGN )
            {
                // ToDo: we could check that an exception is thrown

                continue;
            }
        }

        binary::BinaryOp binop = binary::BinaryOp( op );

        z = applyBinary( x, binop, y );

        BOOST_CHECK( z != T( 0 ) );
    }
}

