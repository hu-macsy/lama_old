/**
 * @file ReductionOpTest.cpp
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
 * @brief Test enum for ReductionOp
 *
 * @author Thomas Brandes
 * @date 30.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/ReductionOp.hpp>
#include <sstream>

using namespace scai;
using namespace common;

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
