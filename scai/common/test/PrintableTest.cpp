/**
 * @file PrintableTest.cpp
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
 * @brief Test functionality of base class Printable
 * @author Thomas Brandes
 * @date 10.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Printable.hpp>

// Printable class uses defaults writeAt

class X1 : public scai::common::Printable
{
};

// Printable class overrides defaults writeAt

class X2 : public scai::common::Printable
{
    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "X2object";
    }
};

BOOST_AUTO_TEST_CASE( PrintableTest )
{
    X1 x1;
    X2 x2;
    std::ostringstream out;
    out << x1;
    BOOST_CHECK_EQUAL( out.str(), typeid( X1 ).name() );
    out.str( "" );
    out.clear();
    out << x2;
    BOOST_CHECK_EQUAL( out.str(), "X2object" );
    scai::common::Printable& p = x2;
    out.str( "" );
    out.clear();
    out << p;
    BOOST_CHECK_EQUAL( out.str(), "X2object" );
    out.str( "" );
    out.clear();
    p.Printable::writeAt( out );
    BOOST_CHECK_EQUAL( out.str(), typeid( X2 ).name() );
}
