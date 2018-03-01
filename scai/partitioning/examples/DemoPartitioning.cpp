/**
 * @file DemoPartitioning.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Demo program for using partitioning
 * @author Thomas Brandes
 * @date 23.08.2017
 */

#include <scai/partitioning.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace partitioning;

typedef DefaultReal ValueType;

int main( int narg, const char* argv[] )
{
    if ( narg < 2 )
    {
        std::cout << "call " << argv[0] << " <matrixFileName>" << std::endl;
    }

    std::string fileName = argv[1];

    auto A = read<CSRSparseMatrix<ValueType>>( fileName );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::cout << *comm << ": Read matrix A = " << A << std::endl;

    PartitioningPtr thePartitioning( Partitioning::create( "METIS" ) );

    if ( !thePartitioning.get() )
    {
        thePartitioning = Partitioning::create( "BLOCK" );
    }

    DistributionPtr dist( thePartitioning->partitionIt( comm, A, 1.0f ) );

    A.redistribute( dist, dist );

    std::cout << *comm << ": Write matrix A = " << A << std::endl;

    A.writeToFile( "matrix%r.mtx" );
}
