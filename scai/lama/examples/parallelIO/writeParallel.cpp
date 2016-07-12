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

using namespace std;

using namespace scai;
using namespace lama;
using namespace dmemo;

static common::scalar::ScalarType getType() 
{
    common::scalar::ScalarType type = common::TypeTraits<double>::stype;
    
    std::string val;
    
    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {   
        scai::common::scalar::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );
        
        if ( env_type == scai::common::scalar::UNKNOWN )
        {   
            std::cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << std::endl;
        }
        
        type = env_type;
    }

    return type;
}

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

    if ( rank == MASTER )
    {
        FileIO::write( owners, fileName );
    }

    // just make sure that no other process starts anything before write is finished

    comm->synchronize();
}

static DistributionPtr readDistribution( const std::string& inFileName )
{
    utilskernel::LArray<IndexType> owners;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getRank() == 0 )
    {
         FileIO::read( owners, inFileName );

         IndexType minId = owners.min();
         IndexType maxId = owners.max();

         // prove:  0 <= minId <= maxId < comm->size() 

         cout << "owner array, size = " << owners.size() << ", min = " << minId << ", max = " << maxId << endl;
    }

    DistributionPtr dist( new GeneralDistribution( owners, comm ) );

    return dist;
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc != 3 )
    {
        cout << "Usage: " << argv[0] << " infile_name outfile_name" << endl;
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

    { 
        std::ostringstream outFileName;

        outFileName << comm->getRank() << "." << comm->getSize() << "." << argv[2];

        matrix.getLocalStorage().writeToFile( outFileName.str() );

        cout << comm << ": written local part of matrix to file " << outFileName << endl;
    }


    printDistribution( matrix.getRowDistribution(), "owners.mtx" );

    DistributionPtr newDist = readDistribution( "owners.mtx" );

    std::cout << "read dist = " << *newDist << std::endl;

    // SCAI_ASSERT_EQUAL( matrix.getRowDistribution(), *newDist, "Error" );
}
