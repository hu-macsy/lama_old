/**
 * @file common/examples/DemoGrid.cpp
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
 * @brief Example of using the Grid classes.
 * @author Thomas Brandes
 * @date 24.05.2017
 */

#include <iostream>
#include <scai/common/Grid.hpp>

/* -----------------------------------------------------------------------------*/

using namespace scai;

int main()
{

    common::Grid1D grid1( 10 );
    common::Grid2D grid2( 2, 5 );
    common::Grid3D grid3( 2, 3, 2 ); 

    std::cout << grid2.linearPos( 1, 4 ) << std::endl;      // prints 9 
    std::cout << grid3.linearPos( 1, 1, 1 ) << std::endl;   // prints 9

    for ( IndexType i = 0; i < grid3.size(); ++i )
    {
        IndexType q[3];
        grid3.gridPos( q, i );   
        std::cout << "Element " << i << " in grid " << grid3 << " is " 
                  << q[0] << ", " << q[1] << ", " << q[2] << std::endl;
    }
}

