/**
 * @file writeParallel.cpp
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
 * @brief Conversion between matrix file formats, uses FileIO factory
 * @author Thomas Brandes
 * @date 19.06.2016
 */

#include <scai/lama/io/PartitionIO.hpp>
#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>
#include <scai/dmemo.hpp>

#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

#include <scai/common/Settings.hpp>

#include <memory>

using namespace std;

using namespace scai;
using namespace lama;
using namespace dmemo;

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    if ( argc < 3 )
    {
        cout << "Usage: " << argv[0] << " infile_name outfile_name" << endl;
        cout << "   file format is chosen by suffix, e.g. frm, mtx, txt, psc, lmf"  << endl;
        cout << "   outfile_name shoud contain %r for partitioned IO"  << endl;
        cout << "   " << endl;
        return -1;
    }

    typedef DefaultReal ValueType;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    CSRSparseMatrix<ValueType> matrix;

    std::string inFileName = argv[1];

    // use supported file format

    matrix.readFromFile( inFileName );

    cout << *comm << ": read CSR matrix : " << matrix << endl;

    auto dist =  dmemo::cyclicDistribution( matrix.getNumRows(), 2, comm );

    matrix.redistribute( dist, matrix.getColDistributionPtr() );

    std::string outFileName = argv[2];

    bool isPartitioned;

    PartitionIO::getPartitionFileName( outFileName, isPartitioned, *comm );

    if ( comm->getSize() > 1 && !isPartitioned )
    {
        std::cout << argv[2] << " is not a partitioned file name" << std::endl;
        return( -1 );
    }
   
    {
        CSRStorage<ValueType>& localPart = matrix.getLocalStorage();

        auto myGlobalRowIndexes = dist->ownedGlobalIndexes();

        IndexType numSet = localPart.setDiagonalFirst( myGlobalRowIndexes );

        bool okay = numSet == localPart.getNumRows();

        okay = comm->all( okay );

        if ( !okay )
        {
            std::cout << *comm << ": there are missing entries on the main diagonal, will stop" << std::endl;
            return( -1 );
        }

        std::cout << *comm << ": write to file " << outFileName 
                  << " my local part: " << localPart << std::endl;

        localPart.writeToFile( outFileName );
    }

    // now read in the partitioned matrix

    auto myPart = read<CSRStorage<ValueType>>( outFileName );

    // make sure that all processors use the same number of columns 

    myPart.resetNumColumns( comm->max( myPart.getNumColumns() ) );

    // get the owned indexes

    hmemo::HArray<IndexType> myOwnedIndexes = myPart.getIA();

    myOwnedIndexes.resize( myPart.getNumRows() );

    utilskernel::HArrayUtils::gather( myOwnedIndexes, myPart.getJA(), myOwnedIndexes, common::BinaryOp::COPY );

    DistributionPtr newDist;

    try
    { 
        newDist = dmemo::generalDistribution( myOwnedIndexes, comm );
    }
    catch( common::Exception& e )
    {
        std::cout << "Failed to build a general distribution: " << e.what() << std::endl;
        std::cout << dist->getCommunicator() << ": myOwnedIndexes = " << hmemo::printIt( myOwnedIndexes ) << std::endl;
        return -1;
    }

    SCAI_ASSERT_ERROR( utilskernel::HArrayUtils::all( myOwnedIndexes, common::CompareOp::EQ, dist->ownedGlobalIndexes() ), "different owners" );

    CSRSparseMatrix<ValueType> matrixNew( dmemo::generalDistribution( std::move( myOwnedIndexes ) ),
                                          std::move( myPart ) );

    SCAI_ASSERT_LT_ERROR( matrix.maxDiffNorm( matrixNew ), 1e-6, "read matrix is different" )
}
