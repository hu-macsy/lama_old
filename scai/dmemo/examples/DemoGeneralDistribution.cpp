/**
 * @file DemoGeneralDistribution.cpp
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
 * @brief Demo program for geneneral distribution
 * @author Thomas Brandes
 * @date 03.02.2017
 */

#include <scai/dmemo.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/tracing.hpp>

using namespace scai;
using namespace dmemo;
using namespace hmemo;
using namespace utilskernel;

using scai::common::BinaryOp;

int main()
{
    SCAI_REGION( "Main.main" )

    SCAI_LOG_THREAD( "Main" )

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType globalSize = 50000;  // number of elements to distribute

    const IndexType np = comm->getSize();  // available processors

    HArray<IndexType> owners( globalSize );

    if ( comm->getRank() == 0 )
    {
        // create random number owners, only host is relevant

        HArrayUtils::setRandom( owners, np - 1 );   // random values from 0 to np - 1

        // HArrayUtils::setScalar( owners, np, BinaryOp::MODULO );

        std::cout << "owners = " << owners << std::endl;

        ReadAccess<IndexType> rOwners( owners );

        for ( IndexType i = 0; i < globalSize; ++i )
        {
            SCAI_ASSERT_VALID_INDEX( rOwners[i], np, "Illegal owner at i = " << i )
        }
    }

    DistributionPtr dist = generalDistributionByOwners( owners, 0, comm );

    std::cout << *comm << ", dist = " << *dist << std::endl;
}
