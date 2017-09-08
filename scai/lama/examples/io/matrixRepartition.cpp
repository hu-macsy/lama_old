/**
 * @file matrixRepartition.cpp
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
 * @brief Repartitioning of a matrix file with n1 partitions into matrix file with n2 partitions
 * @author Thomas Brandes
 * @date 28.10.2016
 */

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>

#include <scai/lama/io/PartitionIO.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/shared_ptr.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/Communicator.hpp>

#include <scai/utilskernel/LArray.hpp>

using namespace std;

using namespace scai;
using namespace hmemo;
using namespace utilskernel;
using namespace lama;

/** Define StoragePtr as shared pointer to any storage */

typedef common::shared_ptr<_MatrixStorage> StoragePtr;

static common::scalar::ScalarType getType()
{
    common::scalar::ScalarType type = common::TypeTraits<double>::stype;

    string val;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        scai::common::scalar::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );

        if ( env_type == scai::common::scalar::UNKNOWN )
        {
            cout << "SCAI_TYPE=" << val << " illegal, is not a scalar type" << endl;
        }

        type = env_type;
    }

    return type;
}

void printHelp( const char* cmd )
{
    cout << "Usage: " << cmd << " [--SCAI_xxx=vvv] infile_name np_in outfile_name np_out" << endl;
    cout << "   infile_name   is input file with the full matrix, must contain %r if np_in > 1"  << endl;
    cout << "   np_in         is the number of partitions used for input matrix"  << endl;
    cout << "   outfile_name  is output file for the partitioned matrices, must contain %r if np_out > 1"  << endl;
    cout << "   np_out        specifies number of partitions for output matrix"  << endl;
    cout << endl;
    cout << "   Note: file format is chosen by suffix, e.g. frm, mtx, txt, psc"  << endl;
    cout << endl;
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
}

/** Read in a matrix from multiple files, each containing a contiguous block of rows.
 *
 *  @param[out] storage is the vertical (row) concatenation of the matrices read from multiple file
 *  @param[in]  inFileName name of the input file, must contain a "%r" to have unique file name for each block
 *  @param[in]  np         number of files among which the matrix is partitioned
 */
void readStorageBlocked( _MatrixStorage& storage, const string& inFileName, const IndexType np )
{
    vector<StoragePtr> storageVector;

    // proof that all input files are available, read them and push the storages

    for ( PartitionId ip = 0; ip < np; ++ ip )
    {
        string inFileNameBlock = inFileName;

        bool isPartitioned;

        PartitionIO::getPartitionFileName( inFileNameBlock, isPartitioned, ip, np );

        SCAI_ASSERT( FileIO::fileExists( inFileNameBlock ),
                     "Input file for block " << ip << " of " << np << " = "
                     << inFileNameBlock << " could not be opened" )

        StoragePtr blockStorage( storage.newMatrixStorage() );   // new temporary storage for each block

        blockStorage->readFromFile( inFileNameBlock );

        storageVector.push_back( blockStorage );
    }

    if ( storageVector.size() == 1 )
    {
        storage.swap( *storageVector[0] );
    }
    else
    {
        IndexType ns = storageVector.size();

        common::scoped_array<const _MatrixStorage*> storages( new const _MatrixStorage*[ ns ] );

        for ( IndexType i = 0; i < ns; ++i )
        {
            storages[i] = storageVector[i].get();
        }

        storage.cat( 0, storages.get(), ns );   // vertical concatenation, dim = 0
    }
}

/** The following method is a special case where the input file contains the full matrix
 *  and where the matrix should be saved into multiple partions
 *
 *  Important: this algorithm does not require memory for the full matrix
 */
void directPartitioning( const string& inFileName, const string& outFileName, const PartitionId np_out )
{
    common::scalar::ScalarType type = getType();

    common::unique_ptr<FileIO>  inputIO ( FileIO::create( FileIO::getSuffix( inFileName ) ) );

    MatrixStorageCreateKeyType key( _MatrixStorage::Format::CSR, type );

    common::unique_ptr<_MatrixStorage> storage( _MatrixStorage::create( key ) );

    IndexType numRows;     // partitioning is done among numbers of rows
    IndexType numColumns;  // dummy here
    IndexType numValues;   // dummy here

    inputIO->readStorageInfo( numRows, numColumns, numValues, inFileName );

    for ( PartitionId ip = 0; ip < np_out; ++ip )
    {
        string outFileNameBlock = outFileName;

        bool isPartitioned;

        PartitionIO::getPartitionFileName( outFileNameBlock, isPartitioned, ip, np_out );

        IndexType lb;   // lower bound for range on a given partition
        IndexType ub;   // upper bound of range for a given partition

        dmemo::BlockDistribution::getLocalRange( lb, ub, numRows, ip, np_out );

        cout << "Matrix stroage block " << ip << " has range " << lb << " - " << ub
             << ", write to file " << outFileNameBlock << endl;

        inputIO->readStorage( *storage, inFileName, lb, ub - lb );

        storage->writeToFile( outFileNameBlock );
    }
}

/** This method saves a matrix into multiple files each file containing a contiguous block of rows.
 *
 *  @param[in] storage is the matrix that is written
 *  @param[in] outFileName name of output file, must contain "%r" as placeholder for block id
 *  @param[in] np_out number of blocks to write
 *
 *  This method is exaclty the same as writing the matrix block distributed with np number of processors.
 */

void writeStorageBlocked( const _MatrixStorage& storage, const string& outFileName, const IndexType np )
{
    // create temporary storage for each block of same type / format

    common::unique_ptr<_MatrixStorage> blockStorage( storage.newMatrixStorage() );

    IndexType numRows = storage.getNumRows();

    for ( PartitionId ip = 0; ip < np; ++ip )
    {
        string outFileNameBlock = outFileName;

        bool isPartitioned;  // here used as dummy

        PartitionIO::getPartitionFileName( outFileNameBlock, isPartitioned, ip, np );

        IndexType lb;   // lower bound for range on a given partition
        IndexType ub;   // upper bound of range for a given partition

        dmemo::BlockDistribution::getLocalRange( lb, ub, numRows, ip, np );

        cout << "Matrix block " << ip << " has range " << lb << " - " << ub
             << ", write to file " << outFileNameBlock << endl;

        if ( np == 1 )
        {
            storage.writeToFile( outFileNameBlock );
        }
        else
        {
            storage.copyBlockTo( *blockStorage, lb, ub - lb );
            blockStorage->writeToFile( outFileNameBlock );
        }
    }
}

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    if ( comm->getSize() > 1 )
    {
        cout << "This program runs only serially, but communicator is " << *comm << endl;
        return -1;
    }

    if ( argc != 5 )
    {
        cout << "Number of arguments is incorrect." << endl << endl;
        printHelp( argv[0] );
        return -1;
    }

    common::scalar::ScalarType type = getType();

    MatrixStorageCreateKeyType key( _MatrixStorage::Format::CSR, type );

    // common::unique_ptr<_MatrixStorage> fullStorage ( _MatrixStorage::create( key ) );
    // common::unique_ptr<_MatrixStorage> blockStorage ( _MatrixStorage::create( key ) );

    // take double as default

    string inFileName  = argv[1];
    string outFileName = argv[3];

    PartitionId np_in  = 0;
    PartitionId np_out = 0;

    istringstream input1( argv[2] );

    input1 >> np_in;

    SCAI_ASSERT( !input1.fail(), "illegal: np_in=" << argv[2] << " in argument list" )

    istringstream input2( argv[4] );

    input2 >> np_out;

    SCAI_ASSERT( !input2.fail(), "illegal: np_out=" << argv[4] << " in argument list" )

    if ( np_in < 1 && np_out > 1 )
    {
        cout << "Partitioning is done block-wise from input file " << inFileName << endl;

        try
        {
            directPartitioning( inFileName, outFileName, np_out );
            return 0;
        }
        catch ( common::Exception& ex )
        {
            cout << "Direct partitioning failed, error: " << ex.what() << endl;
            return -1;
        }
    }

    StoragePtr fullStorage( _MatrixStorage::create( key ) );

    // read in one or all partitions in the memory

    if ( inFileName.find( "%r" ) == string::npos )
    {
        cout << "Read complete storage from single file " << inFileName << endl;

        if ( np_in > 1 )
        {
            cout << "Attention: np_in = " << np_in << " ignored, no %r in filename " << inFileName << endl;
        }

        fullStorage->readFromFile( inFileName );
    }
    else
    {
        // write the storage in np_out block partitions in separate files

        readStorageBlocked( *fullStorage, inFileName, np_in );

        cout << "Storage (merged of " << np_in << " blocks) : " << *fullStorage << endl;
    }

    if ( outFileName.find( "%r" ) == string::npos )
    {
        cout << "Write complete storage to single file " << outFileName << endl;

        if ( np_out > 1 )
        {
            cout << "WARNING: np_out = " << np_out << " ignored, no %r in filename " << outFileName << endl;
        }

        fullStorage->writeToFile( outFileName );
    }
    else
    {
        cout << "Write storage in " << np_out << " blocks to file " << outFileName << endl;

        // write the storage in np_out block partitions in separate files

        writeStorageBlocked( *fullStorage, outFileName, np_out );
    }
}
