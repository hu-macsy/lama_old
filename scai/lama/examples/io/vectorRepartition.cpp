/**
 * @file vectorRepartition.cpp
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
 * @brief Repartitioning of a vector file with n1 partitions into vector file with n2 partitions
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

typedef common::shared_ptr<_HArray> ArrayPtr;

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
    cout << "   infile_name   is input file with the full vector, must contain %r if np_in > 1"  << endl;
    cout << "   np_in         is the number of partitions used for input vector"  << endl;
    cout << "   outfile_name  is output file for the partitioned vector, must contain %r if np_out > 1"  << endl;
    cout << "   np_out        specifies number of partitions for output vector"  << endl;
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

/** The following method is a special case where the input file contains the full vector and
 *  where the vector should be saved into multiple partions 
 * 
 *  Important: this algorithm does not require memory for the full array
 */
void directPartitioning( const string& inFileName, const string& outFileName, const PartitionId np_out )
{
    common::scalar::ScalarType type = getType();

    common::unique_ptr<FileIO>  inputIO ( FileIO::create( FileIO::getSuffix( inFileName ) ) );

    common::unique_ptr<_HArray> array( _HArray::create( type ) );

    IndexType size;    // total size of the array

    inputIO->readArrayInfo( size, inFileName );

    for ( PartitionId ip = 0; ip < np_out; ++ip )
    {
        string outFileNameBlock = outFileName;
  
        bool isPartitioned;

        PartitionIO::getPartitionFileName( outFileNameBlock, isPartitioned, ip, np_out );

        IndexType lb;   // lower bound for range on a given partition
        IndexType ub;   // upper bound of range for a given partition

        dmemo::BlockDistribution::getLocalRange( lb, ub, size, ip, np_out );

        cout << "Vector block " << ip << " has range " << lb << " - " << ub 
             << ", write to file " << outFileNameBlock << endl;

        inputIO->readArrayBlock( *array, inFileName, lb, ub - lb );

        FileIO::write( *array, outFileNameBlock );
    }
}

/** Read in the vector in core form the input file that can be partitioned */

void readVector( _HArray& array, const string& inFileName, const IndexType np_in )
{
    vector<ArrayPtr> storageVector;

    // proof that all input files are available, read them and push the storages

    for ( PartitionId ip = 0; ip < np_in; ++ ip )
    {
        string inFileNameBlock = inFileName;

        bool isPartitioned;

        PartitionIO::getPartitionFileName( inFileNameBlock, isPartitioned, ip, np_in );
   
        SCAI_ASSERT( FileIO::fileExists( inFileNameBlock ), 
                     "Input file for block " << ip << " of " << np_in << " = " 
                     << inFileNameBlock << " could not be opened" )

        ArrayPtr blockVector( _HArray::create( array.getValueType() ) );

        FileIO::read( *blockVector, inFileNameBlock );
  
        storageVector.push_back( blockVector );

        if ( !isPartitioned && np_in > 1 )
        {
            cout << "Input file name does not contain %r for multiple partitions, np_in = " << np_in << " ignored." << endl;
            break;
        }
    }

    if ( storageVector.size() == 1 )
    {
        array.swap( *storageVector[0] );
    }
    else
    {
        // concatenate all input arrays to the result array

        IndexType size = 0;   // determine at first size to get the right size

        for ( size_t i = 0; i < storageVector.size(); ++i )
        {
            size += storageVector[i]->size();
        }

        array.resize( size );

        IndexType offset = 0;

        for ( size_t i = 0; i < storageVector.size(); ++i )
        {
            IndexType localSize = storageVector[i]->size();
            HArrayUtils::setArraySection( array, offset, 1, *storageVector[i], 0, 1, localSize );
            offset += localSize;
        }
    }
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

    vector<ArrayPtr> storageVector;

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

        ArrayPtr blockVector( _HArray::create( type ) );

        FileIO::read( *blockVector, inFileNameBlock );
  
        storageVector.push_back( blockVector );

        if ( !isPartitioned && np_in > 1 )
        {
            cout << "Input file name does not contain %r for multiple partitions, np_in = " << np_in << " ignored." << endl;
            break;
        }
    }

    ArrayPtr fullVector( _HArray::create( type ) );
 
    // read in one or all partitions in the memory

    readVector( *fullVector, inFileName, np_in );

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

            FileIO::write( *fullVector, outFileNameBlock );
            cout << "Write complete vector to file " << outFileNameBlock << endl;
            break;
        }
        else
        {
            if ( !isPartitioned )
            {
                return -1;
            }

            IndexType n = fullVector->size();

            IndexType lb;   // lower bound for range on a given partition
            IndexType ub;   // upper bound of range for a given partition

            dmemo::BlockDistribution::getLocalRange( lb, ub, n, ip, np_out );

            cout << "Vector block " << ip << " has range " << lb << " - " << ub 
                      << ", write to file " << outFileNameBlock << endl;

            ArrayPtr blockVector( _HArray::create( type ) );
            blockVector->resize( ub - lb );
            HArrayUtils::setArraySection( *blockVector, 0, 1, *fullVector, lb, 1, ub - lb );
            FileIO::write( *blockVector, outFileNameBlock );
        }
    }
}
