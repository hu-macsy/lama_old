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

void bench( const IndexType N )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr blockDist( new BlockDistribution( N, comm ) );
    DistributionPtr repDist( new NoDistribution( N ) );

    // generate random numbers

    DenseVector<double> X;
    DenseVector<IndexType> perm;

    float fillRate = 1.0f;

    srand( 131 + comm->getRank() );

    X.setRandom( blockDist, fillRate );

    DenseVector<double> Xrep( X, repDist );   // save unsorted vector

    double tmpTime = Walltime::get();

    bool ascending = false;

    // X.sort( ascending );
  
    X.sort( perm, ascending );

    tmpTime = Walltime::get() - tmpTime;

    std::cout << "Sort time: " << tmpTime << " seconds" << std::endl;

    const utilskernel::LArray<double>& localValues = X.getLocalValues();
    const utilskernel::LArray<IndexType>& permValues = perm.getLocalValues();

    for ( IndexType i = 0; i < X.getDistribution().getLocalSize(); ++i )
    {
        std::cout << "X[local:" << i << "] = " << localValues[i] << std::endl;
        std::cout << "perm[local:" << i << "] = " << permValues[i] << std::endl;
    }

    // check the sorted values

    std::cout << "X.isSorted( " << ascending << " ) = " << X.isSorted( ascending ) << std::endl;

    // check the sorted values

    X.redistribute( repDist );
    perm.redistribute( repDist );
 
    const HArray<double>& repLocalValues = X.getLocalValues();

    SCAI_ASSERT( utilskernel::HArrayUtils::isSorted( repLocalValues, ascending ), "Vector X is not sorted correctly." )

    utilskernel::LArray<double> sortedValues;

    utilskernel::HArrayUtils::gather( sortedValues, Xrep.getLocalValues(), perm.getLocalValues(), utilskernel::binary::COPY );

    std::cout << "diff = " << sortedValues.maxDiffNorm( X.getLocalValues() ) << std::endl;
}

int main()
{
    bench( 20 );
}
