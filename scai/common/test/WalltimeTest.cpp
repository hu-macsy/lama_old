/**
 * @file WalltimeTest.cpp
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
 * @brief Test routines for class Walltime
 * @author Thomas Brandes
 * @date 10.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/common/Walltime.hpp>

BOOST_AUTO_TEST_CASE( WalltimeTest )
{
    int errorPercent = 2;  // should not be less otherwise test might fail from time to time

    using scai::common::Walltime;
    using scai::common::INTEGER_8;

    for ( int k = 0; k < 20; ++k )
    {
        // get some time stamps to avoid possible initialization overhead
        Walltime::timestamp();
    }

    INTEGER_8 i0 = Walltime::timestamp();
    double t0 = Walltime::get();
    Walltime::sleep( 1000 );
    double t1 = Walltime::get();
    INTEGER_8 i1 = Walltime::timestamp();
    // time in seconds
    double time = t1 - t0;
    // should be rather accurate one second, but we give it 2 percent
    BOOST_CHECK_CLOSE( 1.0, time, errorPercent );
    // using timestamp instead of get() should give same result
    BOOST_CHECK_CLOSE( double( i1 - i0 ) / double( Walltime::timerate() ), time, errorPercent );
}
