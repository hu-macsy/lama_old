/**
 * @file PrintableTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Test functionality of base class Printable
 *
 * @author Thomas Brandes
 * @date 10.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Printable.hpp>

#include <unistd.h>

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
    BOOST_CHECK_EQUAL( out.str(), typeid(X1).name() );

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
    BOOST_CHECK_EQUAL( out.str(), typeid(X2).name() );
}
