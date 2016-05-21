/**
 * @file common/examples/DemoLoopMacro.cpp
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
 * @brief Demo of using loop macros to implement nested repeats.
 * @author Thomas Brandes
 * @date 22.05.2016
 */


#include <scai/common/macros/loop.hpp>
#include <iostream>

#define DOIT( x, y ) std::cout << x << y << std::endl;

#define INNER_LOOP( x ) SCAI_COMMON_LOOP_LVL2( x, DOIT, 5, 6, 7 )
#define OUTER_LOOP SCAI_COMMON_LOOP( INNER_LOOP, 1, 2, 3 )

int main()
{
    OUTER_LOOP
}
