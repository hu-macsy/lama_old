/**
 * @file matrixRepartition.cpp
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

#include <scai/utilskernel/LArray.hpp>

using namespace std;

using namespace scai;
using namespace hmemo;
using namespace utilskernel;
using namespace lama;

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

int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

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
 
    if ( np_in < 1 )
    {
        COMMON_THROWEXCEPTION( "Illegal number of partitions for input file = " << np_in )
    }
 
    istringstream input2( argv[4] );

    input2 >> np_out;
 
    if ( np_out < 1 )
    {
        COMMON_THROWEXCEPTION( "Illegal number of partitions for output file = " << np_out )
    }

    typedef common::shared_ptr<_MatrixStorage> StoragePtr;

    vector<StoragePtr> storageVector;

    // proof that all input files are available, read them and push the storages

    for ( PartitionId ip = 0; ip < np_in; ++ ip )
    {
        string inFileNameBlock = inFileName;

        bool isPartitioned;

        PartitionIO::getPartitionFileName( inFileNameBlock, isPartitioned, ip, np_in );
   
        if ( !FileIO::fileExists( inFileNameBlock ) )
        {
            cerr << "Input file for block " << ip << " of " << np_in << " = " << inFileNameBlock << " not available" << endl;
            return -1;
        }

        StoragePtr blockStorage( _MatrixStorage::create( key ) );

        blockStorage->readFromFile( inFileNameBlock );
  
        storageVector.push_back( blockStorage );

        if ( !isPartitioned && np_in > 1 )
        {
            cout << "Input file name does not contain %r for multiple partitions, np_in = " << np_in << " ignored." << endl;
            break;
        }
    }

    StoragePtr fullStorage; 

    if ( storageVector.size() == 1 )
    {
        // avoid concatenation with corresponding copies 

        fullStorage = storageVector[0];
    }
    else
    {
        // concatenate all input storages in a new CSR storage

        fullStorage.reset( _MatrixStorage::create( key ) );
        fullStorage->rowCat( storageVector );
    }

    for ( PartitionId ip = 0; ip < np_out; ++ip )
    {
        string outFileNameBlock = outFileName;
  
        bool isPartitioned;

        PartitionIO::getPartitionFileName( outFileNameBlock, isPartitioned, ip, np_out );

        if ( np_out == 1 || !isPartitioned )
        {
            if ( np_out > 1 )
            {
                cout << "Output file name does not contain %r for multiple partitions, np_out = " << np_out << " ignored" << endl;
            }

            fullStorage->writeToFile( outFileNameBlock );
            cout << "Write complete storage to file " << outFileNameBlock << endl;
            break;
        }
        else
        {
            if ( !isPartitioned )
            {
                return -1;
            }

            IndexType numRows = fullStorage->getNumRows();

            IndexType lb;   // lower bound for range on a given partition
            IndexType ub;   // upper bound of range for a given partition

            dmemo::BlockDistribution::getLocalRange( lb, ub, numRows, ip, np_out );

            cout << "Matrix block " << ip << " has range " << lb << " - " << ub 
                      << ", write to file " << outFileNameBlock << endl;

            StoragePtr blockStorage( _MatrixStorage::create( key ) );
            fullStorage->copyBlockTo( *blockStorage, lb, ub - lb );
            blockStorage->writeToFile( outFileNameBlock );
        }
    }
}
