/**
 * @file elementwiseOpTest.cpp
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
 * @brief Test enum for ElementwiseOp
 * @author Lauretta Schubert
 * @date 06.10.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/ElementwiseOp.hpp>
#include <sstream>

using namespace scai;
using namespace utilskernel;

BOOST_AUTO_TEST_CASE( ElementwiseOpTest )
{
    int count = 0;
    for ( int type = 0; type < elementwise::MAX_ELEMENTWISE_OP; ++type )
    {
        std::ostringstream s;
        s << elementwise::ElementwiseOp( type );
        BOOST_CHECK( s.str().length() > 0 );

        // check if all strings are correct

        if ( type == elementwise::INVERT )
        {
            BOOST_CHECK_EQUAL( s.str(), "INVERT" );
            count++;
        }
        if ( type == elementwise::CONJ )
        {
            BOOST_CHECK_EQUAL( s.str(), "CONJ" );
            count++;
        }
        if ( type == elementwise::EXP )
        {
            BOOST_CHECK_EQUAL( s.str(), "EXP" );
            count++;
        }
        if ( type == elementwise::SQRT )
        {
            BOOST_CHECK_EQUAL( s.str(), "SQRT" );
            count++;
        }
        if ( type == elementwise::SIN )
        {
            BOOST_CHECK_EQUAL( s.str(), "SIN" );
            count++;
        }
        if ( type == elementwise::COS )
        {
            BOOST_CHECK_EQUAL( s.str(), "COS" );
            count++;
        }
        if ( type == elementwise::TAN )
        {
            BOOST_CHECK_EQUAL( s.str(), "TAN" );
            count++;
        }
        if ( type == elementwise::ATAN )
        {
            BOOST_CHECK_EQUAL( s.str(), "ATAN" );
            count++;
        }
        if ( type == elementwise::LOG )
        {
            BOOST_CHECK_EQUAL( s.str(), "LOG" );
            count++;
        }
    }

    // check if all types are tested
    BOOST_CHECK_EQUAL( count, elementwise::MAX_ELEMENTWISE_OP );
}
