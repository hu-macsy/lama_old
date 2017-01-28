/**
 * @file common/examples/DemoFunction.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Example for using function / bind
 * @author Thomas Brandes
 * @date 25.01.2016
 */

#include <scai/common/function.hpp>
#include <scai/common/bind.hpp>

#include <iostream>

int f ( double x, int p )
{
    double f = 1.0;

    for ( int i = 0; i < p; ++i )
    {
        f *= 10.0;
    }

    return static_cast<int>( x * f );
}

using namespace scai::common;

int main()
{
    int ( *foo1 ) ( double, int ) = &f;       // traditional function pointer
    function<int( double, int )> foo2 = &f;    // function wrapper class
    std::cout << "Call foo1( 3.15143, 2 ) = " << foo1( 3.15143, 2 ) << std::endl;
    std::cout << "Call foo2( 3.15143, 3 ) = " << foo2( 3.15143, 3 ) << std::endl;
    function<int( double )> foo2a = bind( f, _1, 5 );
    std::cout << "Call foo2a( 3.15143 ) = " << foo2a( 3.15143 )  << std::endl;
    function<int( int )> foo2b = bind( f, 3.15143, _1 );
    std::cout << "Call foo2b( 1 ) = " << foo2b( 1 )  << std::endl;
    function<int( int, double )> foo3 = bind( f, _2, _1 );
    std::cout << "Call foo3( 2, 3.15143 ) = " << foo3( 2, 3.15143 ) << std::endl;
}
