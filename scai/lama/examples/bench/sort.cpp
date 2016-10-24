/**
 * @file maxnorm.cpp
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
 * @brief Test of maxnorm for all valuetypes
 * @author Eric Schricker
 * @date 21.03.2016
 */

#include <iostream>
#include <iomanip>

#include <scai/lama.hpp>

// Matrix & vector related includes
#include <scai/lama/DenseVector.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/Walltime.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace dmemo;
using namespace std;
using scai::common::Walltime;

int main()
{
    const IndexType N = 20;  // global size of the sorting vector

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr dist( new BlockDistribution( N, comm ) );

    // generate random numbers

    DenseVector<double> X;

    float fillRate = 1.0f;

    srand( 131 + comm->getRank() );

    X.setRandom( dist, fillRate );

    double tmpTime = Walltime::get();

    bool ascending = true;

    X.sort( ascending );

    tmpTime = Walltime::get() - tmpTime;

    std::cout << "Sort time: " << tmpTime << " seconds" << std::endl;

    const utilskernel::LArray<double>& localValues = X.getLocalValues();

    for ( IndexType i = 0; i < X.getDistribution().getLocalSize(); ++i )
    {
        std::cout << "X[local:" << i << "] = " << localValues[i] << std::endl;
    }

    // check the sorted values

    DistributionPtr repDist( new NoDistribution( N ) );

    X.redistribute( repDist );
 
    const HArray<double>& repLocalValues = X.getLocalValues();

    SCAI_ASSERT( utilskernel::HArrayUtils::isSorted( repLocalValues, ascending ), "Vector X is not sorted correctly." )
}
