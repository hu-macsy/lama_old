/**
 * @file Lib/example.cpp
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
 * @brief Simple example executable for testing dllimport/dllexport with
 * @author Thomas Brandes
 * @date 01.05.2013
 */
 		  simple rectangle class.
 * @author Thomas Brandes
 * @date 01.05.2013
 */

#include <iostream>
#include "Rectangle.hpp"

int main( int, char* [] )
{
    double x = 4.0;
    double y = 2.5;

    Rectangle r( x, y );

    std::cout << "Size of rectangle( x = " << x << ", y = " << y << " ) is " << r.getArea() << std::endl;

    return 0;
}
