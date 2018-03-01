/**
 * @file square.cpp
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
 */

#include <scai/lama.hpp>

#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/partitioning/Partitioning.hpp>

#include <memory>

using namespace scai;
using namespace hmemo;
using namespace lama;
using namespace partitioning;

typedef DefaultReal ValueType;

int main( int narg, const char* argv[] )
{
    if ( narg < 3 )
    {
        std::cout << "Please call: " << argv[0] << " <matrixFileName> <numPartitions>" << std::endl;
        return -1;
    }

    std::string fileName = argv[1];
    std::string procs    = argv[2];

    IndexType np = 0;
    HArray<float> pWeights;

    if ( procs.find( "." ) == std::string::npos )
    {
        // assume that procs argument is the number of processors
        np = atoi( argv[2] );
        float weight = static_cast<float>( 1 ) / static_cast<float>( np );
        pWeights.setSameValue( np, weight );
    }
    else
    {
        FileIO::read( pWeights, procs );
        np = pWeights.size();
        // scale the weight:  pWeights /= pWeights.sum();
        WriteAccess<float> rWeight( pWeights );
        float sumWeight = 0;
        for ( IndexType i = 0; i < np; ++i )
        {
            sumWeight += rWeight[ i ];
        }
        for ( IndexType i = 0; i < np; ++i )
        {
            rWeight [ i ] /= sumWeight;
        }
    }

    SCAI_ASSERT_GT_ERROR( np , 0, "Partitioning requires at least one processor" )

    std::string kind = "METIS";
  
    if ( narg > 3 )
    {
        kind = argv[3];
    }
 
    PartitioningPtr partitioning;

    if ( Partitioning::canCreate( kind ) )
    {
        partitioning = Partitioning::create( kind );
    }
    else
    {
        std::cout << "Partitioning kind = " << kind << " not supported." << std::endl;

        for ( auto val : Partitioning::getCreateValues() )
        {
            std::cout << "Supported: " << val << std::endl;
        }

        return EXIT_FAILURE;
    }

    auto csrMatrix = read<CSRSparseMatrix<ValueType>>( fileName );

    HArray<PartitionId> dist;

    partitioning->squarePartitioning( dist, csrMatrix, pWeights );

    std::string distFileName = "dist_" + kind + "_" + std::to_string( np ) + ".txt";

    FileIO::write( dist, distFileName.c_str() );

    std::cout << "Written distribution to " << distFileName << std::endl;
}
