/**
 * @file VectorKindTest.cpp
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
 * @brief Test enum class VectorKind
 * @author Thomas Brandes
 * @date 24.11.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/VectorKind.hpp>
#include <scai/common/test/TestMacros.hpp>
#include <sstream>

using namespace scai;
using namespace lama;

BOOST_AUTO_TEST_CASE( VecorKindTest )
{
    // arithmetic binary operations

    for ( int kind = 0; kind < static_cast<int>( VectorKind::UNDEFINED ); ++kind )
    {
        std::ostringstream s;
        s << VectorKind( kind );
        BOOST_CHECK( s.str().length() > 0 );
    }


    BOOST_CHECK_EQUAL( str2VectorKind( "SPARSE" ), VectorKind::SPARSE );
    BOOST_CHECK_EQUAL( str2VectorKind( "DENSE" ), VectorKind::DENSE );
    BOOST_CHECK_EQUAL( str2VectorKind( "JOINED" ), VectorKind::JOINED );
    BOOST_CHECK_EQUAL( str2VectorKind( "NONSENSE" ), VectorKind::UNDEFINED );
}
