/**
 * @file common/examples/DemoPointer.cpp
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
 * @brief Example of pointer use
 *        /
 * @author Thomas Brandes
 * @date 26.08.2015
 */

#include <iostream>

#include <memory>

/* -----------------------------------------------------------------------------*/

int main()
{
    std::unique_ptr<double> sum( new double );
    std::unique_ptr<double[]> vals ( new double[10] );
    *sum = 0.0;

    for ( int i = 0; i < 10; ++ i )
    {
        vals[i] = i;
    }

    for ( int i = 0; i < 10; ++ i )
    {
        *sum += vals[i];
    }

    std::cout << "Sum = " << *sum << std::endl;
}

