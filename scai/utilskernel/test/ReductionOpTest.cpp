/**
 * @file ReductionOpTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Test enum for ReductionOp
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/ReductionOp.hpp>
#include <sstream>

using namespace scai;
using namespace utilskernel;

BOOST_AUTO_TEST_CASE( ReductionOpTest )
{
    for ( int type = reduction::COPY; type <= reduction::ABS_MAX + 1; ++type )
    {
        std::ostringstream s;
        s << reduction::ReductionOp( type );
        BOOST_CHECK( s.str().length() > 0 );
        if ( type == reduction::COPY )
        {
            BOOST_CHECK_EQUAL( s.str(), "COPY" );
        }
    }
}
