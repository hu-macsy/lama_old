/**
 * @file ScalarTypeTest.cpp
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
 * @brief Test enum for ScalarType
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/mepr/ScalarTypeHelper.hpp>
#include <scai/common/SCAITypes.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_CASE( ScalarTypeTest )
{
    for ( int type = 0; type <= static_cast<int>( ScalarType::UNKNOWN ); ++type )
    {
        ScalarType stype = ScalarType( type );
        std::ostringstream s;
        s << stype;
        BOOST_CHECK( s.str().length() > 0 );
        BOOST_CHECK_EQUAL( stype, str2ScalarType( s.str().c_str() ) );

#ifdef SCAI_COMPLEX_SUPPORTED

        // This test will fail as TypeTraits for Complex are not available

        size_t pos = s.str().find( "Complex" );

        if ( isComplex( stype ) )
        {
            BOOST_CHECK( pos != std::string::npos );
        }
        else
        {
            BOOST_CHECK( pos == std::string::npos );
        }
#endif
    }

    BOOST_CHECK( typeSize( ScalarType::PATTERN ) == 0 );

    BOOST_CHECK( ! isNumeric( ScalarType::INT ) );

    BOOST_CHECK( ! isNumeric( ScalarType::INDEX_TYPE ) );

    BOOST_CHECK_EQUAL( sizeof( IndexType ), typeSize( ScalarType::INDEX_TYPE ) );

    BOOST_CHECK_THROW (
    {
        typeSize( ScalarType::INTERNAL );
    }, common::Exception );

    bool contains1 = mepr::ScalarTypeHelper<SCAI_TYPELIST( char, int )>::contains( ScalarType::INT );
    BOOST_CHECK( contains1 );

    bool contains2 = mepr::ScalarTypeHelper<SCAI_TYPELIST( char, int )>::contains( ScalarType::FLOAT );
    BOOST_CHECK( !contains2 );
}

BOOST_AUTO_TEST_CASE( precisionTest )
{
    for ( int type = 0; type < static_cast<int>( ScalarType::UNKNOWN ); ++type )
    {
        ScalarType stype = ScalarType( type );

        if ( stype == ScalarType::INTERNAL )
        {
            // should throw an exception to make sure that call is replaced with correct type

            BOOST_CHECK_THROW(
            {
                precision( stype );
            }, common::Exception );

            continue;
        }

        int n = precision( stype );

        if ( isNumeric( stype ) )
        {
            BOOST_CHECK( n > 0 );
        }
        else
        {
            BOOST_CHECK_EQUAL( 0, n );
        }
    }
}


