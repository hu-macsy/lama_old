/**
 * @file redistribute.cpp
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
 * @brief ToDo: Missing description in ./partitioning/examples/redistribute.cpp
 * @author Thomas.Brandes@scai.fraunhofer.de 2018-03-01
 * @date 16.03.2015
 */

#include <scai/lama.hpp>

#include <scai/dmemo/RedistributePlan.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai;
using namespace lama;

int main( int, char** )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType N = 10;

    dmemo::DistributionPtr sourceDistribution( new dmemo::BlockDistribution( N, comm ) );
   
    typedef DefaultReal ValueType;

    DenseVector<ValueType> v = linearDenseVector( sourceDistribution, ValueType( 1 ), ValueType( 2 ) );

    std::cout << "v = " << v << std::endl;

    hmemo::HArray<PartitionId> newLocalOwners;

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x = v[i];
        std::cout << "v[ " << i << " ] = " << x << std::endl;
    }

    { 
        IndexType nLocal = sourceDistribution->getLocalSize();
        IndexType npart = comm->getSize();

        hmemo::WriteOnlyAccess<PartitionId> wMapping( newLocalOwners, nLocal );

        for ( IndexType i = 0; i < nLocal; ++i )
        {
            IndexType globalI = sourceDistribution->local2Global( i );
            IndexType globalChunk = globalI / 3;
            wMapping[i] = globalChunk % npart;
        }
    }

    auto redist = dmemo::redistributePlanByNewOwners( newLocalOwners, sourceDistribution );

    v.redistribute( redist );

    std::cout << "v = " << v << std::endl;

    for ( IndexType i = 0; i < N; ++i )
    {
        ValueType x = v[i];
        std::cout << "v[ " << i << " ] = " << x << std::endl;
    }
}
