/**
 * @file Lib/example.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief ToDo: Missing description in ./Lib/example.cpp
 * @author Thomas Brandes
 * @date 01.05.2013
 */
// Simple example program that prints size of pointer (32-bit or 64-bit)

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
