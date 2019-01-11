/**
 * @file BinaryOpTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Test enum for BinaryOp
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/BinaryOp.hpp>
#include <scai/common/CompareOp.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_SUITE( BinaryOpTest )

BOOST_AUTO_TEST_CASE( BinaryOpTest )
{
    // arithmetic binary operations

    for ( int type = 0; type < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++type )
    {
        std::ostringstream s;
        s << BinaryOp( type );
        BOOST_CHECK( s.str().length() > 0 );

        if ( BinaryOp( type ) == BinaryOp::COPY )
        {
            BOOST_CHECK_EQUAL( s.str(), "COPY" );
        }
    }

    // comparison binary operations

    for ( int type = 0; type < static_cast<int>( CompareOp::MAX_COMPARE_OP ); ++type )
    {
        std::ostringstream s;
        s << CompareOp( type );
        BOOST_CHECK( s.str().length() > 0 );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( ApplyBinOpTest, T, scai_array_test_types )
{
    // we do not check here the results but make sure that all ops are supported

    for ( int op = 0; op < static_cast<int>( BinaryOp::MAX_BINARY_OP ); ++op )
    {
        T x = 5;
        T y = 2;
        T z = 0;

        if ( common::TypeTraits<T>::stype == common::TypeTraits<IndexType>::stype )
        {
            if ( BinaryOp( op ) == BinaryOp::POW || BinaryOp( op ) == BinaryOp::COPY_SIGN )
            {
                // ToDo: we could check that an exception is thrown

                continue;
            }
        }

        BinaryOp binop = BinaryOp( op );

        z = applyBinary( x, binop, y );

        BOOST_CHECK( z != T( 0 ) );
    }

    // test compare operations

    BOOST_CHECK( compare<T>( T( 2 ), CompareOp::LT, T( 3 ) ) );
    BOOST_CHECK( !compare<T>( T( 3 ), CompareOp::LT, T( 3 ) ) );
    BOOST_CHECK( compare<T>( T( 3 ), CompareOp::LE, T( 3 ) ) );
    BOOST_CHECK( compare<T>( T( 3 ), CompareOp::GE, T( 3 ) ) );
    BOOST_CHECK( !compare<T>( T( 3 ), CompareOp::GT, T( 3 ) ) );
    BOOST_CHECK( compare<T>( T( 3 ), CompareOp::GT, T( 2 ) ) );
}

BOOST_AUTO_TEST_SUITE_END()
