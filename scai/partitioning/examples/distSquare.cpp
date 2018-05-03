/**
 * @file matrix.cpp
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
 * @brief Construction of matrix and vector with move operations
 * @author Thomas Brandes
 * @date 2018-01-02 
 */

#include <scai/lama.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/partitioning/ParMetisPartitioning.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Walltime.hpp>

#include <iostream>
#include <cstdlib>

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace dmemo;
using partitioning::PartitioningPtr;
using partitioning::Partitioning;

typedef DefaultReal ValueType;

int main ( int argc, char* argv[] )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // Important: First region should be entered after communicator has been build

    SCAI_REGION( "user.main" )

    if ( argc < 2 )
    {   
        std::cerr << "No input file specified" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Read a sparse matrix from the file that has been specified by command line argument

    auto matrix = read<CSRSparseMatrix<ValueType>>( argv[1] );

    std::cout << "Read matrix: " << matrix << std::endl;
  
    SCAI_ASSERT_EQ_ERROR( matrix.getNumRows(), matrix.getNumColumns(), "matrix not square" )

    // make initial distribution, can be arbitrary

    auto dist = std::make_shared<CyclicDistribution>( matrix.getNumRows(), 2, comm );
 
    {
        SCAI_REGION( "user.distribute" )
        matrix.redistribute( dist, dist );
    }

    std::cout << "Initial matrix: " << matrix << std::endl;

    // now compute new owners with Partioning

    PartitioningPtr myPartioning = Partitioning::create( "PARMETIS" );

    float weight = 1.0f;

    HArray<IndexType> newLocalOwners;

    {
        SCAI_REGION( "user.partitioning" )
        myPartioning->squarePartitioning( newLocalOwners, matrix, weight );
    }

    std::string filename = "newOwners_" + std::to_string( comm->getRank() ) + ".txt";

    FileIO::write( newLocalOwners, filename );

    // Now build a Redistributor that contains permutation and communicaton pattern

    double time = common::Walltime::get();

    SCAI_REGION_START( "user.buildRedistributor" )
    dmemo::Redistributor redist( newLocalOwners, dist );
    SCAI_REGION_END( "user.buildRedistributor" )

    time = common::Walltime::get() - time;

    DistributionPtr newDist = redist.getTargetDistributionPtr();


    std::cout << "new target distribution = " << *newDist << std::endl;
    std::cout << "time for building redist = " << time << std::endl;

    // Multiply with old matrix

    std::cout << "Matrix multiplication with original matrix: " << matrix << std::endl;

    auto x  = linearDenseVector( dist, ValueType( 1 ), ValueType( 1.0 ) );
    auto y1 = eval<DenseVector<ValueType>>( matrix * x );

    std::cout << "Initial vector x = " << x << std::endl;

    time = common::Walltime::get();

    x.redistribute( redist );

    time = common::Walltime::get() - time;

    std::cout << "Redistributed vector x = " << x << std::endl;
    std::cout << "Time for redistribution x = " << time << std::endl;

    time = common::Walltime::get();

    {
        SCAI_REGION( "user.redistMatrix" )
        matrix.redistribute( redist, newDist );
    }

    time = common::Walltime::get() - time;

    std::cout << "time for redistribute matrix: " << time << std::endl;

    std::cout << "Matrix-vector multiplication with redistributed matrix: " << matrix << std::endl;

    auto y2 = eval<DenseVector<ValueType>>( matrix * x );

    y2.writeToFile( "y2.txt" );

    y1.redistribute( newDist );

    ValueType diff = y1.maxDiffNorm( y2 );
 
    std::cout << "max diff = " << diff << std::endl;
}
