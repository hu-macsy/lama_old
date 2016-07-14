/**
 * @file writeParallel.cpp
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
 * @brief Conversion between matrix file formats, uses FileIO factory
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>
#include <scai/dmemo.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/unique_ptr.hpp>

#include "utility.hpp" 

using namespace std;

using namespace scai;
using namespace lama;
using namespace dmemo;

static void printDistribution( const Distribution& distribution, const std::string& fileName )
{
    using namespace hmemo;

    const PartitionId MASTER = 0;

    CommunicatorPtr comm = distribution.getCommunicatorPtr();
 
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    if ( size == 1 )
    {
        return;   // do not print a NoDistribution
    }

    HArray<IndexType> indexes;

    if ( rank == MASTER )
    {
        // we need the owners only on the host processor
        // indexes = 0, 1, 2, ..., globalSize - 1

        utilskernel::HArrayUtils::setOrder( indexes, distribution.getGlobalSize() );
    }

    HArray<IndexType> owners;

    // Note: only master process asks for owners, other processes have 0 indexes

    distribution.computeOwners( owners, indexes );

    std::cout << *comm << ", owner computation finished, owners = " << owners << std::endl;

    if ( rank == MASTER )
    {
        std::cout << *comm << ", MASTER, write distribution to " << fileName << std::endl;

        FileIO::write( owners, fileName );
    }

    // just make sure that no other process starts anything before write is finished

    comm->synchronize();
}

static void printPDistribution( const Distribution& distribution, const std::string& fileName )
{
    // each processor writes a file with its global indexes

    using namespace hmemo;

    CommunicatorPtr comm = distribution.getCommunicatorPtr();
 
    const IndexType nLocal = distribution.getLocalSize();
    const IndexType nGlobal = distribution.getGlobalSize();

    HArray<IndexType> myGlobalIndexes;

    {
        WriteOnlyAccess<IndexType> wGlobalIndexes( myGlobalIndexes, nLocal );

        IndexType k = 0;

        for ( IndexType i = 0; i < nGlobal; ++i )
        {
            if ( distribution.isLocal( i ) )
            {
                wGlobalIndexes[k++] = i;
            }
        }
        
        SCAI_ASSERT_EQ_ERROR( k, nLocal, "serious local mismatch" );
    }

    FileIO::write( myGlobalIndexes, fileName );
}

static DistributionPtr readDistribution( const std::string& inFileName )
{
    utilskernel::LArray<IndexType> owners;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getRank() == 0 )
    {
        std::cout << *comm << ", MASTER, read distribution from " << inFileName << std::endl;

        if ( FileIO::fileExists( inFileName ) )
        {
            FileIO::read( owners, inFileName );

            IndexType minId = owners.min();
            IndexType maxId = owners.max();

             // prove:  0 <= minId <= maxId < comm->size() 
    
            cout << "owner array, size = " << owners.size() << ", min = " << minId << ", max = " << maxId << endl;
        }
    }

    IndexType ownersSum = comm->sum( owners.size() );

    DistributionPtr dist( new GeneralDistribution( owners, comm ) );

    if ( ownersSum > 0 )
    {
        dist.reset( new GeneralDistribution( owners, comm ) );
    }
    else
    {
        dist.reset( new BlockDistribution( 0, comm ) );
    }

    return dist;
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        cout << "Usage: " << argv[0] << " infile_name outfile_name distfile_name" << endl;
        cout << "   file format is chosen by suffix, e.g. frm, mtx, txt, psc"  << endl;
        cout << "   --SCAI_TYPE=<data_type> is data type of input file and used for internal representation" << endl;
        cout << "   --SCAI_IO_BINARY=0|1 to force formatted or binary output file" << endl;
        cout << "   --SCAI_IO_TYPE_DATA=<data_type> is data type used for file output" << endl;
        cout << "   " << endl;
        cout << "   Supported types: ";
        vector<common::scalar::ScalarType> dataTypes;
        hmemo::_HArray::getCreateValues( dataTypes );
        for ( size_t i = 0; i < dataTypes.size(); ++i )
        { 
            cout << dataTypes[i] << " ";
        }
        cout << endl;
        return -1;
    }

    // take double as default 

    common::scalar::ScalarType type = getType();

    // oops, no factory for storage, only for matrix

    common::unique_ptr<Matrix> matrixPtr( Matrix::getMatrix( Matrix::CSR, type ) );

    Matrix& matrix = *matrixPtr;

    std::string inFileName = argv[1];

    // use supported file format

    matrix.readFromFile( inFileName );

    cout << "read CSR matrix : " << matrix << endl;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr dist( new dmemo::CyclicDistribution( matrix.getNumRows(), 2, comm ) );

    matrix.redistribute( dist, matrix.getColDistributionPtr() );

    std::string outFileName = argv[2];

    bool writePartitions;

    getPartitionFileName( outFileName, writePartitions, *comm );

    if ( !writePartitions )
    {
        // write it in one single file 

        matrix.writeToFile( outFileName );

        cout << comm << ": written matrix to file " << outFileName << endl;
    }
    else
    { 
        matrix.getLocalStorage().writeToFile( outFileName );

        cout << *comm << ": written local part of matrix to file " << outFileName << endl;
    }

    if ( argc > 3 )
    {
        std::string distFileName = argv[3];

        getPartitionFileName( distFileName, writePartitions, *comm );

        if ( writePartitions )
        {
            printPDistribution( matrix.getRowDistribution(), distFileName );
        }
        else
        {
            printDistribution( matrix.getRowDistribution(), distFileName );
            DistributionPtr newDist = readDistribution( distFileName );
            SCAI_ASSERT_EQ_ERROR( newDist->getGlobalSize(), matrix.getNumRows(), "mismatch" )
        }

    }
}
