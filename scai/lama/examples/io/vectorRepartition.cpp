/**
 * @file vectorRepartition.cpp
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
 * @brief Repartitioning of a vector file with n1 partitions into vector file with n2 partitions
 * @author Thomas Brandes
 * @date 28.10.2016
 */

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama.hpp>

#include <scai/lama/io/PartitionIO.hpp>

#include <scai/common/Settings.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <memory>

using namespace std;

using namespace scai;
using namespace hmemo;
using namespace utilskernel;
using namespace lama;

static common::ScalarType getType()
{
    common::ScalarType type = common::TypeTraits<double>::stype;

    string val;

    if ( scai::common::Settings::getEnvironment( val, "SCAI_TYPE" ) )
    {
        scai::common::ScalarType env_type = scai::common::str2ScalarType( val.c_str() );

        if ( env_type == scai::common::ScalarType::UNKNOWN )
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
    vector<common::ScalarType> dataTypes;
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
    common::ScalarType type = getType();

    std::unique_ptr<FileIO>  inputIO ( FileIO::create( FileIO::getSuffix( inFileName ) ) );

    std::unique_ptr<_HArray> array( _HArray::create( type ) );

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

        inputIO->readArray( *array, inFileName, lb, ub - lb );

        FileIO::write( *array, outFileNameBlock );
    }
}

/** Read in the vector in core form the input file that can be partitioned */

void readArrayBlocked( _HArray& array, const string& inFileName, const IndexType np_in )
{
    typedef std::shared_ptr<_HArray> ArrayPtr;    // use shared pointer in vector

    vector<ArrayPtr> blockVector;

    // proof that all input files are available, read them and push the storages

    for ( PartitionId ip = 0; ip < np_in; ++ ip )
    {
        string inFileNameBlock = inFileName;

        bool isPartitioned;

        PartitionIO::getPartitionFileName( inFileNameBlock, isPartitioned, ip, np_in );

        SCAI_ASSERT( FileIO::fileExists( inFileNameBlock ),
                     "Input file for block " << ip << " of " << np_in << " = "
                     << inFileNameBlock << " could not be opened" )

        ArrayPtr blockArray( _HArray::create( array.getValueType() ) );

        FileIO::read( *blockArray, inFileNameBlock );

        blockVector.push_back( blockArray );
    }

    if ( blockVector.size() == 1 )
    {
        array.swap( *blockVector[0] );
    }
    else
    {
        // concatenate all input arrays to the result array

        IndexType size = 0;   // determine at first size to get the right size

        for ( size_t i = 0; i < blockVector.size(); ++i )
        {
            size += blockVector[i]->size();
        }

        array.resize( size );

        IndexType offset = 0;

        for ( size_t i = 0; i < blockVector.size(); ++i )
        {
            IndexType localSize = blockVector[i]->size();
            HArrayUtils::_setArraySection( array, offset, 1, *blockVector[i], 0, 1, localSize );
            offset += localSize;
        }
    }
}

/** This method saves a vector into multiple files each file containing a contiguous block of values
 *
 *  @param[in] array is the vector that is written
 *  @param[in] outFileName name of output file, must contain "%r" as placeholder for block id
 *  @param[in] np_out number of blocks to write
 *
 *  This method is exaclty the same as writing the vector block distributed with np number of processors.
 */

void writeArrayBlocked( const _HArray& array, const string& outFileName, const IndexType np )
{
    // create temporary array for each block of same type

    std::unique_ptr<_HArray> blockArray( _HArray::create( array.getValueType() ) );

    IndexType size = array.size();

    for ( PartitionId ip = 0; ip < np; ++ip )
    {
        string outFileNameBlock = outFileName;

        bool isPartitioned;  // here used as dummy

        PartitionIO::getPartitionFileName( outFileNameBlock, isPartitioned, ip, np );

        IndexType lb;   // lower bound for range on a given partition
        IndexType ub;   // upper bound of range for a given partition

        dmemo::BlockDistribution::getLocalRange( lb, ub, size, ip, np );

        cout << "Matrix block " << ip << " has range " << lb << " - " << ub
             << ", write to file " << outFileNameBlock << endl;

        if ( np == 1 )
        {
            FileIO::write( array, outFileNameBlock );
        }
        else
        {
            blockArray->clear();  // invalidate content to avoid unnecessary mem transfer
            blockArray->resize( ub - lb );
            HArrayUtils::_setArraySection( *blockArray, 0, 1, array, lb, 1, ub - lb );
            FileIO::write( *blockArray, outFileNameBlock );
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

    common::ScalarType type = getType();

    MatrixStorageCreateKeyType key( Format::CSR, type );

    // std::unique_ptr<_MatrixStorage> fullStorage ( _MatrixStorage::create( key ) );
    // std::unique_ptr<_MatrixStorage> blockStorage ( _MatrixStorage::create( key ) );

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

    std::unique_ptr<_HArray> fullVector( _HArray::create( type ) );

    // read in one or all partitions in the memory

    if ( inFileName.find( "%r" ) == string::npos )
    {
        cout << "Read complete storage from single file " << inFileName << endl;

        if ( np_in > 1 )
        {
            cout << "Attention: np_in = " << np_in << " ignored, no %r in filename " << inFileName << endl;
        }

        FileIO::read( *fullVector, inFileName );
    }
    else
    {
        // read in one or all partitions in the memory

        readArrayBlocked( *fullVector, inFileName, np_in );

        cout << "Array (merged of " << np_in << " blocks) : " << *fullVector << endl;
    }

    if ( outFileName.find( "%r" ) == string::npos )
    {
        cout << "Write complete vector to single file " << outFileName << endl;

        if ( np_out > 1 )
        {
            cout << "WARNING: np_out = " << np_out << " ignored, no %r in filename " << outFileName << endl;
        }

        FileIO::write( *fullVector, outFileName );
    }
    else
    {
        cout << "Write vector in " << np_out << " blocks to file " << outFileName << endl;

        // write the array in np_out block partitions in separate files

        writeArrayBlocked( *fullVector, outFileName, np_out );
    }
}
