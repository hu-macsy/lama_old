/**
 * @file ScalarTypeTest.cpp
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
 * @brief Test enum for ScalarType
 *
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/ScalarType.hpp>
#include <sstream>

using namespace scai;
using namespace common;

BOOST_AUTO_TEST_CASE( ScalarTypeTest )
{
    for ( int type = scalar::INDEX_TYPE; type <= scalar::UNKNOWN; ++type )
    {
        scalar::ScalarType stype = scalar::ScalarType( type );

        std::ostringstream s;
        s << stype;

        BOOST_CHECK( s.str().length() > 0 );

        BOOST_CHECK_EQUAL( stype, str2ScalarType( s.str().c_str() ) );

        size_t pos = s.str().find( "Complex" );

        if ( isComplex( stype ) )
        {
           BOOST_CHECK( pos != std::string::npos );
        }
        else
        {
           BOOST_CHECK( pos == std::string::npos );
        }
    }
}
