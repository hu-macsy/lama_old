/**
 * @file common/test/UnaryOpTest.cpp
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
 * @brief Test enum for UnaryOp
 * @author Lauretta Schubert
 * @date 06.10.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/UnaryOp.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_CASE( UnaryOpTest )
{
    int count = 0;

    for ( int type = 0; type < static_cast<int>( UnaryOp::MAX_UNARY_OP ); ++type )
    {
        UnaryOp op = UnaryOp( type );

        std::ostringstream s;
        s << op;
        BOOST_CHECK( s.str().length() > 0 );

        // check if all strings are correct

        if ( op == UnaryOp::COPY )
        {
            BOOST_CHECK_EQUAL( s.str(), "COPY" );
            count++;
        }

        if ( op == UnaryOp::CONJ )
        {
            BOOST_CHECK_EQUAL( s.str(), "CONJ" );
            count++;
        }

        if ( op == UnaryOp::SQR )
        {
            BOOST_CHECK_EQUAL( s.str(), "SQR" );
            count++;
        }

        if ( op == UnaryOp::MINUS )
        {
            BOOST_CHECK_EQUAL( s.str(), "MINUS" );
            count++;
        }

        if ( op == UnaryOp::ABS )
        {
            BOOST_CHECK_EQUAL( s.str(), "ABS" );
            count++;
        }

        if ( op == UnaryOp::ASUM )
        {
            BOOST_CHECK_EQUAL( s.str(), "ASUM" );
            count++;
        }

        if ( op == UnaryOp::EXP )
        {
            BOOST_CHECK_EQUAL( s.str(), "EXP" );
            count++;
        }

        if ( op == UnaryOp::SQRT )
        {
            BOOST_CHECK_EQUAL( s.str(), "SQRT" );
            count++;
        }

        if ( op == UnaryOp::SIN )
        {
            BOOST_CHECK_EQUAL( s.str(), "SIN" );
            count++;
        }

        if ( op == UnaryOp::COS )
        {
            BOOST_CHECK_EQUAL( s.str(), "COS" );
            count++;
        }

        if ( op == UnaryOp::TAN )
        {
            BOOST_CHECK_EQUAL( s.str(), "TAN" );
            count++;
        }

        if ( op == UnaryOp::ATAN )
        {
            BOOST_CHECK_EQUAL( s.str(), "ATAN" );
            count++;
        }

        if ( op == UnaryOp::LOG )
        {
            BOOST_CHECK_EQUAL( s.str(), "LOG" );
            count++;
        }

        if ( op == UnaryOp::FLOOR )
        {
            BOOST_CHECK_EQUAL( s.str(), "FLOOR" );
            count++;
        }

        if ( op == UnaryOp::CEIL )
        {
            BOOST_CHECK_EQUAL( s.str(), "CEIL" );
            count++;
        }

        if ( op == UnaryOp::SIGN )
        {
            BOOST_CHECK_EQUAL( s.str(), "SIGN" );
            count++;
        }

        if ( op == UnaryOp::RECIPROCAL )
        {
            BOOST_CHECK_EQUAL( s.str(), "RECIPROCAL" );
            count++;
        }

    }

    // check if all types are tested
    BOOST_CHECK_EQUAL( UnaryOp( count ), UnaryOp::MAX_UNARY_OP );
}
