/**
 * @file TestDistributions.hpp
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
 * @brief Class that provides a vector of different distributions
 * @author Thomas Brandes
 * @date 25.07.2016
 */

#pragma once

#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/SingleDistribution.hpp>

namespace scai
{

namespace dmemo
{

/* --------------------------------------------------------------------- */

/** Class that is a vector of different distribution pointers.
 *
 */
class TestDistributions : public std::vector<DistributionPtr>
{
public:

    /** Constructor of different test distributions.
     *
     *  @param[in] globalSize is same global size for all distributions
     */
    TestDistributions( const IndexType globalSize )
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        // Take all available distributions from the factory

        std::vector<std::string> values;

        Distribution::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

            BOOST_CHECK_EQUAL( dist->getKind(), values[i] );

            if ( values[i] == "METIS" )
            {
                // METIS prints a lot of warnings with a single processor, so skip it

                if ( comm->getSize() == 1 )
                {
                    continue;
                }

                // METIS prints warnings if number of elements is less than global size

                if ( static_cast<IndexType>( comm->getSize() ) > globalSize )
                {
                    continue;
                }
            }

            push_back( dist );
        }

        // Create a random general distribution, must be same on all processors

        hmemo::HArray<PartitionId> owners;

        {
            PartitionId owner = 315;
            PartitionId numPartitions = comm->getSize();

            hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, globalSize );

            for ( IndexType i = 0; i < globalSize; ++i )
            {
                owner = owner * 119 % 185;
                wOwners[i] = owner % numPartitions;
            }
        }

        PartitionId root = 0;

        push_back( generalDistributionByOwners( owners, root, comm ) );

        // Create a general block distribution with different weights on each processor

        float weight = static_cast<float>( comm->getRank() + 1 );

        push_back( genBlockDistributionByWeight( globalSize, weight, comm ) );

        // Create a single distributon, not on first processor
        //  1 -> 0, 2 ->1, 3 -> 1, 4 -> 2

        PartitionId owner = comm->getSize() / 2;

        push_back( DistributionPtr( new SingleDistribution( globalSize, comm, owner ) ) );
    }

private:

    TestDistributions();
};

/* --------------------------------------------------------------------- */

}

}



