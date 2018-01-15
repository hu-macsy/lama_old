/**
 * @file ContextTypeTest.cpp
 * Copyright (c) 2009-2017
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
 * @brief Test enum for different context types.
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/ContextType.hpp>
#include <scai/common/AccessKind.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_CASE( ContextTypeTest )
{
    for ( int type = 0; type <= static_cast<int>( ContextType::MaxContext ); ++type )
    {
        std::ostringstream s;
        s << ContextType( type );
        BOOST_CHECK( s.str().length() > 0 );

        if ( ContextType( type ) == ContextType::Host )
        {
            BOOST_CHECK_EQUAL( s.str(), "Host" );
        }
    }
}

BOOST_AUTO_TEST_CASE( AccessKindTest )
{
    for ( int kind = 0; kind <= static_cast<int>( AccessKind::MaxAccessKind ); ++kind )
    {
        std::ostringstream s;
        s << AccessKind( kind );
        BOOST_CHECK( s.str().length() > 0 );

        if ( AccessKind( kind ) == AccessKind::Read )
        {
            // output should contain at least an R for read and no W at all
            BOOST_CHECK( s.str().find( "R" ) != std::string::npos );
            BOOST_CHECK( s.str().find( "W" ) == std::string::npos );
        }
    }
}

