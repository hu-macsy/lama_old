/**
 * @file StorageIO.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of static IO routines for matrix storage
 * @author Thomas Brandes
 * @date 27.07.2012
 */

#include <scai/common/OpenMP.hpp>
#include <scai/common/ScalarType.hpp>

// hpp
#include <scai/lama/StorageIO.hpp>

// local library
#include <scai/lama/io/FileStream.hpp>
//#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/io/mmio.hpp>
#include <scai/lama/io/IOUtils.hpp>
#include <scai/lama/mepr/IOWrapper.hpp>


// internal scai libraries
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;
using common::unique_ptr;
using common::scoped_array;
using sparsekernel::OpenMPCSRUtils;

namespace lama
{

#define VERSION_ID 22

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( _StorageIO::logger, "StorageIO" )

/* -------------------------------------------------------------------------- */

const int _StorageIO::mIversion = 4;

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromSAMGFile(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    IndexType& numColumns,
    const std::string& filename )
{
    SCAI_REGION( "StorageIO.readCSRFromSAMGFile" )

    // start with reading the header
    FileStream inFile( filename, std::ios::in );

    int iversion;
    char fileType = '!';
    IndexType numValues, numRows, id, size, rank;

    numColumns = 0;
    numValues = 0;

    inFile >> fileType >> iversion;

    if( fileType != 'f' && fileType != 'b' )
    {
        COMMON_THROWEXCEPTION( "Invalid file type" )
    }

    if( iversion != mIversion )
    {
        COMMON_THROWEXCEPTION( "Invalid file version: " << iversion << ", should be " << mIversion )
    }

    inFile >> numValues;
    inFile >> numRows;
    // AMG is always square?
    numColumns = numRows;

    inFile >> id;
    inFile >> size;
    inFile >> rank;

    inFile.close(); // explicitly, otherwise done by destructor

    // now read *.amg file in correct mode
    std::ios::openmode flags = std::ios::in;

    if( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    if( filename.size() < 4 )
    {
        COMMON_THROWEXCEPTION( "Invalid filename, can't load *.amg file" )
    }
    std::string filenameData = filename.substr( 0, filename.size() - 4 ) + ".amg";
    inFile.open( filenameData, flags );

    //TODO: allow different type to be read?
    inFile.read( csrIA, numRows + 1, -1, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrJA, numValues, -1, common::TypeTraits<IndexType>::stype, '\n' );
    inFile.read( csrValues, numValues, 0, common::TypeTraits<ValueType>::stype, '\n' );

    inFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToSAMGFile(
    const PartitionId size,
    const PartitionId rank,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& filename,
    const common::scalar::ScalarType iaType,
    const common::scalar::ScalarType jaType,
    const common::scalar::ScalarType valuesType,
    const bool writeBinary )
{
    SCAI_REGION( "StorageIO.writeCSRToBinaryFile " )

    char fileType;
    if( writeBinary )
    {
        fileType = 'b';
    }
    else
    {
        fileType = 'f';
    }

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    FileStream outFile( filename + ".frm", std::ios::out | std::ios::trunc );
    outFile << fileType << " \t" << mIversion << "\n";
    outFile << "\t\t" << numValues << "\t" << numRows << "\t" << VERSION_ID << "\t" << size << "\t" << rank;
    outFile.close();

    SCAI_LOG_INFO( logger, "writeCSRToBinaryFile ( " << filename << ".amg" << ")" << ", #rows = " << csrIA.size()-1
                           << ", #values = " << csrJA.size() )

    std::ios::openmode flags = std::ios::out | std::ios::trunc;
    if( writeBinary )
    {
       flags |= std::ios::binary;
    }

    outFile.open( filename + ".amg", flags );
    outFile.write<IndexType>( csrIA, 1, jaType, '\n' );
    outFile.write<IndexType>( csrJA, 1, iaType, '\n' );
    outFile.write<ValueType>( csrValues, 0, valuesType, '\n' );
    outFile.close();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToMMFile(
    const HArray<IndexType>& csrIA,
    const IndexType numColumns,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName,
    const common::scalar::ScalarType& dataType )
{
    SCAI_REGION( "StorageIO.writeCSRToMMFile" )

    // TODO: different output type possible?

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    FileStream outFile( fileName + ".mtx", std::ios::out | std::ios::trunc );

    if( dataType != common::scalar::INTERNAL )
    {
        writeMMHeader( false, numRows, numColumns, numValues, outFile, dataType );
    }
    else
    {
        writeMMHeader( false, numRows, numColumns, numValues, outFile, common::TypeTraits<ValueType>::stype );
    }


    // output code runs only for host context
    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( csrIA, host );
    ReadAccess<IndexType> ja( csrJA, host );
    ReadAccess<ValueType> data( csrValues, host );

    for( IndexType ii = 0; ii < numRows; ++ii )
    {
        for( IndexType jj = ia[ii]; jj < ia[ii + 1]; ++jj )
        {
            outFile << ii + 1 << " " << ja[jj] + 1;

            if( dataType != common::scalar::PATTERN )
            {
                outFile << " " << data[jj];
            }

            outFile << std::endl;
        }
    }

    outFile.close();
}

/* -------------------------------------------------------------------------- */

// Template struct describes one Value in a Sparse Matrix
// used for the mtx Reader
template<typename ValueType>
struct MatrixValue
{
    IndexType i;
    IndexType j;
    ValueType v;

    MatrixValue( IndexType row, IndexType col, ValueType val )
                    : i( row ), j( col ), v( val )
    {
    }
};

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromMMFile(
    HArray<IndexType>& csrIA,
    IndexType& numColumns,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName )
{
    bool isSymmetric, isPattern;
    IndexType numRows, numValues, numValuesFile;

    FileStream inFile( fileName, std::ios::in );

    readMMHeader( numRows, numColumns, numValuesFile, isPattern, isSymmetric, inFile );

    ContextPtr host = Context::getHostPtr();

    WriteOnlyAccess<IndexType> ia( csrIA, host, numRows + 1 );

    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<IndexType> > setVal;
    setVal.getSupportedContext( host );

    // initialize ia;
    setVal[host]( ia, numRows+1, 0, utilskernel::reduction::COPY );

    std::vector<MatrixValue<ValueType> > values;

    // Set the right size of the Vector
    if( isSymmetric )
    {
        values.reserve( numValuesFile * 2 - numRows );
    }
    else
    {
        values.reserve( numValuesFile );
    }

    //Create Input Vector
    MatrixValue<ValueType> val( 0, 0, 0 );

    SCAI_LOG_DEBUG( logger, "beginning read in" )

    std::string line;

    numValues = numValuesFile;
    // TODO: there could be comment or blank lines in the file (allowed by specification?), we should over-read them!
    for( int l = 0; l < numValuesFile && !inFile.eof(); ++l )
    {
    	std::getline(inFile, line);
    	std::istringstream reader(line);

        reader >> val.i;
        reader >> val.j;

        if( !isPattern )
        {
            reader >> val.v;
        }
        else
        {
            val.v = ValueType( 1.0 );
        }

        ++ia[val.i-1];

        // if the matrix is symmetric, the value appears in row 'column' again.
        if( isSymmetric )
        {
            if( val.j != val.i )
            {
                ++ia[val.j-1];
                ++numValues;

                MatrixValue<ValueType> val_symmetric( val.j - 1, val.i - 1, val.v );
                values.push_back( val_symmetric );
            }
        }

        // TODO: do this in parallel after reading is finsihed?
        val.i--;
        val.j--;
        values.push_back( val );
    }

    if( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline(inFile, line);
    if( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << fileName << "': invalid file, contains to many elements." )
    }

    inFile.close();
    // Create csrJA, csrValues
    WriteOnlyAccess<IndexType> ja( csrJA, numValues );
    WriteOnlyAccess<ValueType> data( csrValues, numValues );

    //create absolute Values of ia
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::scan<IndexType> > scan;
    scan.getSupportedContext( host );

    // convert sizes to offset array
    scan[host]( ia, numRows + 1 );

    //initialize ia and data
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::setVal<ValueType> > setValValues;
    setValValues.getSupportedContext( host );

    setVal[host]( ja, numValues, -1, utilskernel::reduction::COPY );
    setValValues[host]( data, numValues, 0, utilskernel::reduction::COPY );


    for( IndexType elem = 0; elem < numValues; elem++ )
    {
        MatrixValue<ValueType> value = values[elem];
        IndexType offset = ia[value.i];
        IndexType pos = 0;

        while( ja[offset + pos] != -1 )
        {
            pos++;
        }

        SCAI_LOG_TRACE( logger, "added row " << value.i << ", offset = " << offset + pos << ", j = " << value.j )
        ja[offset + pos] = value.j;
        data[offset + pos] = value.v;
    }


    static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;
    sortRowElements.getSupportedContext( host );

    // Note: we do care if the matrix has really all diagonal elements available
    sortRowElements[host]( ja, data, ia, numRows, true );

    SCAI_LOG_INFO( logger, "construct matrix " << numRows << " x " << numColumns << " from CSR arrays, # non-zeros = " << numValues )
}

/* -------------------------------------------------------------------------- */

bool _StorageIO::fileExists( const std::string& fileName )
{
    // open the file for reading -> it exists

    std::ifstream file( fileName.c_str(), std::ios::in );
    return file.good();
}

/* -------------------------------------------------------------------------- */

bool _StorageIO::hasSuffix( const std::string& fileName, const std::string& suffix)
{
    return fileName.size() >= suffix.size() &&
           fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) == 0;
}

/* -------------------------------------------------------------------------- */

void _StorageIO::writeMMHeader(
		const bool& vector,
		const IndexType& numRows,
		const IndexType& numColumns,
		const IndexType& numValues,
		FileStream& outFile,
		const common::scalar::ScalarType& dataType )
{
    outFile << "%%matrixmarket ";

    if( vector )
    {
        outFile << "vector array ";
    }
    else
    {
        outFile << "matrix coordinate ";
    }

    switch( dataType )
    {
        case common::scalar::DOUBLE:
        case common::scalar::FLOAT:
            outFile << "real ";
            break;
        case common::scalar::COMPLEX:
        case common::scalar::DOUBLE_COMPLEX:
        case common::scalar::LONG_DOUBLE_COMPLEX:
            outFile << "complex ";
            break;
        case common::scalar::INDEX_TYPE:
            outFile << "integer ";
            break;
        case common::scalar::PATTERN:
            outFile << "pattern ";
            break;
        default:
            COMMON_THROWEXCEPTION( "_StorageIO::writeMMHeader: " "unknown datatype." << dataType )

    }

    // TODO: Add support for symmetric matrices
    // currently we can only write non-symmetric
    outFile << "general" << std::endl;

    if( vector )
    {
        outFile << numRows << " " << numColumns << std::endl;
    }
    else
    {
        outFile << numRows << " " << numColumns << " " << numValues << std::endl;
    }
}

/* -------------------------------------------------------------------------- */

void _StorageIO::readMMHeader(
		IndexType& numRows,
		IndexType& numColumns,
		IndexType& numValues,
		bool& isPattern,
		bool& isSymmetric,
		FileStream& inFile )
{
    std::string buffer;

    // read %%MatrixMarket
    std::getline(inFile, buffer, ' ' );
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    if( buffer != "%%matrixmarket" )
    {
        COMMON_THROWEXCEPTION( "Given file is no valid matrix market file, expected file to begin with %%MatrixMarket" )
    }

    // read object type
    std::getline(inFile, buffer, ' ' );
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    // check if object type is valid in general
    if( buffer != "matrix" && buffer != "vector" )
    {
        COMMON_THROWEXCEPTION( "Object type in the given matrix market file is invalid, should be matrix or vector" )
    }
    bool isVector = false;
    if( buffer == "vector" )
    {
        isVector = true;
    }


    // read file type
    std::getline(inFile, buffer, ' ' );
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    // checkif file type is valid in general
    if( buffer != "coordinate" && buffer != "array" )
    {
        COMMON_THROWEXCEPTION( "Format type in the given matrix market file is invalid, should be coordinate or array" )
    }
//    if( buffer == "array" )
//    {
//        // TODO: array => dense data, do we need to check this later?
//        COMMON_THROWEXCEPTION( "File type 'array' (dense data) is currently not supported!" )
//    }

    // read data type
    std::getline(inFile, buffer, ' ' );
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    if( buffer != "real" && buffer != "integer" && buffer != "complex" && buffer != "pattern" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be real, integer, complex or pattern" )
    }
    // TODO: allow to return other value types as well => check if the valid type is used
    if( buffer == "pattern" )
    {
        isPattern = true;
    }
    else
    {
        isPattern = false;
    }

    // read symmetry
    std::getline(inFile, buffer, '\n' );
    std::transform(buffer.begin(), buffer.end(), buffer.begin(), ::tolower);
    if( buffer != "general" && buffer != "symmetric" && buffer != "skew-symmetric" && buffer != "hermitian" )
    {
        COMMON_THROWEXCEPTION( "Data type in the given matrix market file is invalid, should be general, symmetric, skew-symmetric or hermitian" )
    }
    if( buffer == "general" )
    {
        isSymmetric = false;
    }
    else
    {
        if( buffer == "symmetric" )
        {
            isSymmetric = true;
        }
        else
        {
            // TODO: add support!
            COMMON_THROWEXCEPTION( "Symmetry options 'skew-symmetric' and 'hermitian' are currently not supported!" )
        }
    }

    do
    {
        std::getline(inFile, buffer, '\n' );
    } while( buffer.at( 0 ) == '%' );

    std::stringstream bufferSS( buffer );
    bufferSS >> numRows;
    bufferSS >> numColumns;
    // TODO: vector correct here? should it be dense vs sparse?
    if( !isVector )
    {
        bufferSS >> numValues;
    }
    else
    {
        numValues = numRows * numColumns;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StorageIO<ValueType>::writeCSRToFile(
    const PartitionId size,
    const PartitionId rank,
    const HArray<IndexType>& csrIA,
    const IndexType numColumns,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const std::string& fileName,
    const File::FileType& fileType,
    const common::scalar::ScalarType& valuesType,
    const common::scalar::ScalarType& iaType,
    const common::scalar::ScalarType& jaType,
    const bool writeBinary /* = false */ )
{
    SCAI_REGION( "StorageIO.writeCSRToFile " )

    std::string fileBaseName;

    if( hasSuffix( fileName, ".frm" ) )
    {
        fileBaseName = fileName.substr( 0, fileName.size() - 4 ).c_str();
    }
    else
    {
        fileBaseName = fileName.c_str();
    }

    // for multiple parititions we generate an own id for each partition

    if( size > 1 )
    {
        char rankstr[10];
        sprintf( rankstr, ".%d", rank );
        fileBaseName += rankstr;
    }

    switch( fileType )
    {
        case File::SAMG:
            writeCSRToSAMGFile( size, rank, csrIA, csrJA, csrValues, fileBaseName, iaType, jaType, valuesType, writeBinary );
            break;
        case File::MATRIX_MARKET:
            writeCSRToMMFile( csrIA, numColumns, csrJA, csrValues, fileBaseName, valuesType );
            return;
        default:
            COMMON_THROWEXCEPTION( "Unknown file type definition." )
    }

}

template<typename ValueType>
void StorageIO<ValueType>::readCSRFromFile(
    HArray<IndexType>& csrIA,
    IndexType& numColumns,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const std::string& fileName )
{
    SCAI_LOG_INFO( logger, "read CSR matrix data from file: '" << fileName << "'." )

    SCAI_REGION( "StorageIO.readCSRFromFile" )

    std::string suffix;
    std::string baseFileName = fileName;
    std::string amgFileName;
    File::FileType fileType;

    if( fileName.size() >= 4 )
    {
        suffix = fileName.substr( fileName.size() - 4, 4 );
    }

    if( suffix == ".frm" )
    {
        fileType = File::SAMG;
    }

    if( suffix == ".mtx" )
    {
    	fileType = File::MATRIX_MARKET;
    }


    switch( fileType )
    {
        case File::SAMG:
            readCSRFromSAMGFile( csrIA, csrJA, csrValues, numColumns, fileName );
            break;

        case File::MATRIX_MARKET:
            readCSRFromMMFile( csrIA, numColumns, csrJA, csrValues, fileName );
            break;

        default:
            COMMON_THROWEXCEPTION( "Read storage file: unknown file type = " << fileType )
    }
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromFile( HArray<ValueType>& data,
                                              IndexType& numColumns,
                                              const std::string& filename )
{
    File::FileType fileType;
    std::string suffix;
    std::string used_filename;

    if( filename.size() >= 4 )
    {
        suffix = filename.substr( filename.size() - 4, 4 );
    }

    if( suffix == ".frm" )
    {
        fileType = File::SAMG;
        used_filename = filename.substr( 0, filename.size() - 4 ) + ".frv";
    }

    if( suffix == ".frv" )
    {
        fileType = File::SAMG;
        used_filename = filename;
    }

    if( suffix == ".mtx" )
    {
        fileType = File::MATRIX_MARKET;
        used_filename = filename;
    }

    switch( fileType )
    {
        case File::SAMG:
            readDenseFromSAMGFile( data, numColumns, used_filename );
            break;
        case File::MATRIX_MARKET:
            readDenseFromMMFile( data, numColumns, used_filename );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unknown File Type." );
    }

}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToFile( const HArray<ValueType>& data,
                                             const IndexType& numColumns,
                                             const std::string& filename,
                                             const File::FileType fileType,
                                             const common::scalar::ScalarType dataType,
                                             const bool writeBinary /* = false */ )
{
    switch( fileType )
    {
        case File::SAMG:
            writeDenseToSAMGFile( data, numColumns, filename, dataType, writeBinary );
            break;

        case File::MATRIX_MARKET:
            writeDenseToMMFile( data, numColumns, filename, dataType, writeBinary );
            break;

        default:
            COMMON_THROWEXCEPTION( "Unknown file type definition." );
    }
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromSAMGFile( HArray<ValueType>& data,
                                                  IndexType& numColumns,
                                                  const std::string& filename )
{
    SCAI_LOG_INFO( logger, "read DenseVector from file '" << filename << "'." )

    char fileType;
    int dataTypeSize;
    IndexType numRows;

    // start with reading the *.frv header file
    FileStream inFile( filename, std::ios::in );

    inFile >> fileType;
    inFile >> numRows;
    inFile >> dataTypeSize;
    inFile.close();

    // *.frv files can only store vectors
    numColumns = 1;

    if ( fileType != 'b' && fileType != 'f' )
    {
        COMMON_THROWEXCEPTION( "Invalid header file." )
    }

    // now read *.vec file in correct mode
    std::ios::openmode flags = std::ios::in;

    if( fileType == 'b' )
    {
        flags |= std::ios::binary;
    }

    if( filename.size() < 4 )
    {
        COMMON_THROWEXCEPTION( "Invalid filename, can't load *.vec file" )
    }
    std::string filenameData = filename.substr( 0, filename.size() - 4 ) + ".vec";
    inFile.open( filenameData, flags );

    common::scalar::ScalarType dataType;
    switch(dataTypeSize){
        // special cases for handling IndexType, int and long, as these are not properly supported yet
        case 4:
            dataType = common::scalar::FLOAT;
            break;
        case 8:
            dataType = common::scalar::DOUBLE;
            break;

        default:
			dataType = common::scalar::UNKNOWN;
            SCAI_LOG_ERROR( logger, "Encountered invalid type size " << dataTypeSize )
    }


    inFile.read( data, numRows, 0, dataType, '\n' );

    inFile.close();
}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToSAMGFile( const HArray<ValueType>& data,
                                                 const IndexType& numColumns,
                                                 const std::string& filename,
                                                 const common::scalar::ScalarType dataType,
                                                 const bool writeBinary /* = false */ )
{
    SCAI_ASSERT_ERROR( numColumns == 1, "SAMG format can only store dense vectors" )

    // start by writing the *.frv header file
    char fileType;

    if( writeBinary )
    {
        fileType = 'b';
    }
    else
    {
        fileType = 'f';
    }

    // TODO: maybe provide this as function at a central place?

    int typeSize = common::mepr::ScalarTypeHelper<SCAI_ARITHMETIC_ARRAY_HOST_LIST>::sizeOf( dataType );

    if( typeSize == 0 )
    {
        switch( dataType )
        {
            case common::scalar::INTERNAL:
                typeSize = sizeof(ValueType);
                break;
            default:
                SCAI_LOG_ERROR( logger, "Encountered invalid scalar type " << dataType )
                break;
        }
    }

    FileStream outFile( filename + ".frv", std::ios::out );

    outFile << fileType << std::endl;
    outFile << data.size() << std::endl;
    outFile << typeSize;
    outFile.close();

    // write data into *.vec file
    std::ios::openmode flags = std::ios::out | std::ios::trunc;
    if( writeBinary )
    {
        flags |= std::ios::binary;
    }

    outFile.open( filename + ".vec", flags );

    outFile.write<ValueType>( data, 0, dataType, '\n' );

    outFile.close();
}

template<typename ValueType>
void StorageIO<ValueType>::readDenseFromMMFile( HArray<ValueType>& data,
                                                IndexType& numColumns,
                                                const std::string& filename )
{
    bool isSymmetric, isPattern;
    IndexType numRows, numValues, i;
    ValueType val;
    std::string line;

    FileStream inFile( filename, std::ios::in );

    readMMHeader( numRows, numColumns, numValues, isPattern, isSymmetric, inFile );

    WriteOnlyAccess<ValueType> vector( data, numValues );
    ValueType* vPtr = vector.get();

    for( int l = 0; l < numValues && !inFile.eof(); ++l )
    {
        std::getline( inFile, line );
        std::istringstream reader( line );

        if( isPattern )
        {
            reader >> i;
            val = 1.0;
            i--;
        }
        else
        {
            reader >> val;
            i = l;
        }

        vPtr[i] = val;

    }

    if( inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << filename << "': reached end of file, before having read all data." )
    }

    // check if there is more data in the file tht should not be there
    std::getline(inFile, line);
    if( !inFile.eof() )
    {
        COMMON_THROWEXCEPTION( "'" << filename << "': invalid file, contains to many elements." )
    }

    inFile.close();
    SCAI_LOG_INFO( logger, "construct vector " << numRows )
}

template<typename ValueType>
void StorageIO<ValueType>::writeDenseToMMFile( const HArray<ValueType>& data,
                                               const IndexType& numColumns,
                                               const std::string& filename,
                                               const common::scalar::ScalarType dataType,
                                               const bool writeBinary /* = false */ )
{
    SCAI_ASSERT_ERROR( writeBinary == false, "Matrix market format can not be written binary" );

    FileStream outFile( filename + ".mtx", std::ios::out | std::ios::trunc );

    if( dataType != common::scalar::INTERNAL )
    {
        writeMMHeader( true, data.size(), numColumns, -1, outFile, dataType );
    }
    else
    {
        writeMMHeader( true, data.size(), numColumns, -1, outFile, common::TypeTraits<ValueType>::stype );
    }

    // output code runs only for host context
    ContextPtr host = Context::getHostPtr();

    ReadAccess<ValueType> dataRead( data, host );

    for( IndexType i = 0; i < data.size(); ++i )
    {
        outFile << dataRead[i] << std::endl;
    }

    outFile.close();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( StorageIO, SCAI_ARITHMETIC_HOST )

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
