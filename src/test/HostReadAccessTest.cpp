/**
 * @file HostReadAccessTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class HostReadAccessTest.
 * @author: Alexander Büchel
 * @date 16.07.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <test/TestMacros.hpp>

#include <lama/HostReadAccess.hpp>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HostReadAccessTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.HostReadAccessTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    const IndexType n = 10;

    LAMAArray<double> lamaArray( n, 1.0 );

    HostReadAccess<double> lamaArrayReadAccess( lamaArray );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( 1.0, lamaArrayReadAccess[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    const IndexType n = 10;

    LAMAArray<double> lamaArray( n, 1.0 );

    HostReadAccess<double> lamaArrayRAccess( lamaArray );

    LAMA_WRITEAT_TEST( lamaArrayRAccess );

    lamaArrayRAccess.release();

    LAMA_WRITEAT_TEST( lamaArrayRAccess );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
