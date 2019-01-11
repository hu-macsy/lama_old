/**
 * @file common/test/MatrixOpTest.cpp
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
 * @brief Test enum for MatrixOp
 * @author Thomas Brandes
 * @date 01.03.2018
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/MatrixOp.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_CASE( MatrixOpTest )
{
    for ( int type = 0; type < static_cast<int>( MatrixOp::MAX_MATRIX_OP ); ++type )
    {
        MatrixOp op = MatrixOp( type );

        std::ostringstream s;
        s << op;
        BOOST_CHECK( s.str().length() > 0 );
    }
}

BOOST_AUTO_TEST_CASE( combineTest )
{
    for ( int i1 = 0; i1 < static_cast<int>( MatrixOp::MAX_MATRIX_OP ); ++i1 )
    {
        MatrixOp op1 = MatrixOp( i1 );
   
        BOOST_CHECK_EQUAL( combine( op1, op1 ), MatrixOp::NORMAL );

        for ( int i2 = 0; i2 < static_cast<int>( MatrixOp::MAX_MATRIX_OP ); ++i2 )
        {
            MatrixOp op2 = MatrixOp( i2 );

            BOOST_CHECK_EQUAL( combine( op1, op2 ), combine( op2, op1 ) );
        }
    }
}
