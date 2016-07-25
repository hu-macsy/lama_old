/**
 * @file TestDistributions.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Class that provides a vector of different distributions
 * @author Thomas Brandes
 * @date 25.07.2016
 */

#pragma once

#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/utilskernel.hpp>

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

        std::vector<std::string> values;

        Distribution::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

            BOOST_CHECK_EQUAL( dist->getKind(), values[i] );

            push_back( dist );
        } 

        utilskernel::LArray<PartitionId> owners;

        {
            PartitionId owner = 315;
            PartitionId nPartitions = comm->getSize();

            hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, globalSize );

            for ( IndexType i = 0; i < globalSize; ++i )
            {
                owner = owner * 119 % 185;
                wOwners[i] = owner % nPartitions;
            }
        }

        push_back( DistributionPtr( new GeneralDistribution( owners, comm ) ) );

        float weight = static_cast<float>( comm->getRank() + 1 );

        push_back( DistributionPtr( new GenBlockDistribution( globalSize, weight, comm ) ) );
    }

private:

    TestDistributions();
};

/* --------------------------------------------------------------------- */

}

}



