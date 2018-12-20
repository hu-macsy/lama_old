/**
 * @file partitioning/examples/redistribute1.cpp
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
 * @brief Example of using a redistribute plan for redistribution of a vector.
 * @author Thomas Brandes
 * @date 10.05.2016
 */

#include <scai/lama.hpp>

#include <scai/dmemo/RedistributePlan.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

using namespace scai;
using namespace lama;

int main( int, char** )
{
    const IndexType N = 10;

    auto sourceDistribution = dmemo::blockDistribution( N );
    auto targetDistribution = dmemo::cyclicDistribution( N, 3 );
   
    typedef DefaultReal ValueType;

    auto v = linearDenseVector<ValueType>( sourceDistribution, 1, 0.2 );

    std::cout << "v = " << v << std::endl;

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x = v[i];
        std::cout << "v[ " << i << " ] = " << x << std::endl;
    }

    dmemo::RedistributePlan redist( targetDistribution, sourceDistribution );

    v.redistribute( redist );

    std::cout << "v = " << v << std::endl;

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x = v[i];
        std::cout << "v[ " << i << " ] = " << x << std::endl;
    }
}
