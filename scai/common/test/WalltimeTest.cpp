/**
 * @file WalltimeTest.cpp
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
 * @brief Test routines for class Walltime
 *
 * @author Thomas Brandes
 * @date 10.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Walltime.hpp>

#include <unistd.h>

BOOST_AUTO_TEST_CASE( WalltimeTest )
{

    using scai::common::Walltime;
    using scai::common::INTEGER_8;

    INTEGER_8 i0 = Walltime::timestamp();
    double t0 = Walltime::get();

    sleep( 1 );
    double t1 = Walltime::get();
    INTEGER_8 i1 = Walltime::timestamp();

    // time in seconds

    double time = t1 - t0;

    // should be rather accurate one second

    BOOST_CHECK_CLOSE( 1.0, time, 1 );

    // using timestamp instead of get() should give same result

    BOOST_CHECK_CLOSE( double( i1 - i0 ) / double( Walltime::timerate() ), time, 1 );
}
