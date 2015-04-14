/**
 * @file TracingTest.hpp
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
 * @brief Test for multiple LAMA contexts and transfering context array data
 * @author Alexander Büchel, Thomas Brandes
 * @date 03.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <lama/tracing.hpp>

#ifdef _WIN32
#include <Windows.h>
#endif //WIN32
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TracingTest )

/* --------------------------------------------------------------------- */

void work( int n )
{
#ifdef _WIN32
    Sleep( n );
#else
    sleep( n );
#endif
}
;

void bar()
{
    LAMA_REGION_START( "bar" )
    work( 1 );
    LAMA_REGION_END( "bar" )
}
;

void foo( bool call )
{
    LAMA_REGION( "foo" )
    work( 1 );

    if ( call )
    {
        bar();
        work( 1 );
    }
}
;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( regionTest )
{
    foo( false );
    bar();
    foo( true );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
