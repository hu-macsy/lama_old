/**
 * @file simple.cpp
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
 * @brief simple.cpp
 * @author The LAMA development team
 * @date 17.05.2013
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace lama;

/** Take default real type for this example. */
typedef DefaultReal ValueType;

int main()
{
    //
    // Create a DenseVector of size 8 with value 1.1 in each position
    //
    const IndexType size = 8;

    const DenseVector<ValueType> v( size, 1.1 );
    //
    // Compute the L1 norm of the vector and print it
    //
    auto norm = v.l1Norm();
    std::cout << "L1 norm of v = " << norm << std::endl;
    //
    //  That's it.
    //
    std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;
    return EXIT_SUCCESS;
}
