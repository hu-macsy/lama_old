/**
 * @file common/examples/TimePrecision.cpp
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
 * @brief Check resolution of walltime
 * @author Thomas Brandes
 * @date 10.12.2015
 */

#include <scai/common/Walltime.hpp>

#include <iostream>

using scai::common::Walltime;

int main()
{
    const int NITER = 4;

    for ( int i = 0; i < NITER; ++i )
    {
        long counter = 0;     // count changes of output value of Walltime.:get
        double t = Walltime::get();
        bool stop = false;   // set it to true after one second
        double tc = t;

        while ( !stop )
        {
            double t1 = Walltime::get();

            if ( t1 > tc )
            {
                // value has changed
                counter++;
                tc = t1;
            }

            stop = ( t1 - t ) >= 1.0;
        }

        std::cout << "Resolution: at least " << counter << " ticks per seconds" << std::endl;
    }
}
