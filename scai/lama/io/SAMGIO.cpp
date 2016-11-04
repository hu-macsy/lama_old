/**
 * @file SAMGIO.cpp
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
 * @brief Implementation of IO methods for SAMG format
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#include "SAMGIO.hpp"

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/LArray.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/ScalarType.hpp>

#include <cstring>

/** SAMG file suffixes
 * 
 *  Note: static variables can cause problems as values are already needed during static initialization. 
 */

#define SAMG_MAT_HEADER_SUFFIX ".frm"
#define SAMG_MAT_DATA_SUFFIX   ".amg"
#define SAMG_VEC_HEADER_SUFFIX ".frv"
#define SAMG_VEC_DATA_SUFFIX   ".vec"

#define SAMG_VERSION_ID 22
#define SAMG_IVERSION   4

namespace scai
{

using namespace hmemo;

namespace lama
{


std::string SAMGIO::getVectorFileSuffix() const
{
    return SAMG_VEC_HEADER_SUFFIX;
}
    
std::string SAMGIO::getMatrixFileSuffix() const
{
    return SAMG_MAT_HEADER_SUFFIX;
}

/* --------------------------------------------------------------------------------- */
/*    Implementation of Factory methods                                              */
/* --------------------------------------------------------------------------------- */

FileIO* SAMGIO::create()
{
    return new SAMGIO();
}

std::string SAMGIO::createValue()
{
    return SAMG_MAT_HEADER_SUFFIX;
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeAt( std::ostream& stream ) const
{
    stream << "SAMGIO ( ";
    stream << "suffix = " << getMatrixFileSuffix() << "|" << getVectorFileSuffix() << ", ";
    writeMode( stream );
    stream << " )";
}

/* --------------------------------------------------------------------------------- */

bool SAMGIO::isSupportedMode( const FileMode ) const
{
    // all file modes are supported

    return true;
}

/* --------------------------------------------------------------------------------- */

/** Help routine to get data file name by header file name 
 *
 *  @param[in] headerFileName is the file name of header file
 *  @return    name of the data file
 *
 *  Note: returns same name if no distinction between header and data 
 */
static std::string getDataFileName( const std::string& headerFileName )
{
    std::string result = headerFileName;

    if ( FileIO::hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX) )
    {
        size_t len = strlen( SAMG_MAT_HEADER_SUFFIX );
        result.replace( result.length() - len, len, SAMG_MAT_DATA_SUFFIX );
    }
    else if ( FileIO::hasSuffix( headerFileName, SAMG_VEC_HEADER_SUFFIX ) )
    {
        size_t len = strlen( SAMG_VEC_HEADER_SUFFIX );
        result.replace( result.length() - len, len, SAMG_VEC_DATA_SUFFIX );
    }

    return result;   // same name if no distinction between header and data
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( SAMGIO::logger, "FileIO.SAMGIO" )

/* --------------------------------------------------------------------------------- */

SAMGIO::SAMGIO() 
{
    if ( common::Settings::getEnvironment( mAppendMode, "SCAI_IO_APPEND" ) )
    {
        if ( mAppendMode )
        {
            SCAI_LOG_WARN( logger, "SAMG format does not support append mode" ) 
        }
    }
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeVectorHeader( 
    const IndexType n,
    const IndexType typeSize,
    const bool binary, 
    const std::string& fileName )
{
    char fileType = binary ? 'b' : 'f';

    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    outFile << fileType << std::endl;
    outFile << n << std::endl;
    outFile << typeSize;

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeArrayImpl(
    const hmemo::HArray<ValueType>& array,
    const std::string& fileName )
{
    // needed for header file: type size is size of data type used in output 

    int typeSize = sizeof( ValueType );

    if ( mScalarTypeData != common::scalar::INTERNAL )
    {
        typeSize = common::typeSize( mScalarTypeData );
    }

    bool binary = mFileMode != FORMATTED; 

    writeVectorHeader( array.size(), typeSize, binary, fileName );

    // write data into SAMG vector data file

    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( binary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    IOStream outFile( dataFileName, flags );
 
    if ( binary )
    {
        outFile.writeBinary( array, mScalarTypeData );
    }
    else
    {
        int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

        outFile.writeFormatted( array, precData );
    }

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readVectorHeader( IndexType& n, IndexType& typeSize, bool& binary, const std::string& fileName )
{
    char fileType = ' ';

    typeSize = 0;
    n        = 0;

    IOStream inFile( fileName, std::ios::in );

    inFile >> fileType;
    inFile >> n;
    inFile >> typeSize;

    if ( inFile.fail() )
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << fileName
                               << ", could not read '[f|b] <n> <typeSize>" )
    }

    inFile.close();

    if ( fileType == 'b' )
    {
        binary = true;
    }
    else if ( fileType == 'f' )
    {
        binary = false;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << fileName
                               << ", fileType = " << fileType << " illegal" )
    }
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readArrayInfo( IndexType& size, const std::string& fileName )
{
    IndexType dataTypeSize;   // dummy variable needed for readVectorHeader
    bool binary;              // dummy variable needed for readVectorHeader

    readVectorHeader( size, dataTypeSize, binary, fileName );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readArrayImpl( hmemo::HArray<ValueType>& array, const std::string& fileName, const IndexType first, const IndexType n )
{
    IndexType dataTypeSize = 0;
    IndexType size = 0;
    bool binary = false;

    // start with reading the *.frv header file
 
    readVectorHeader( size, dataTypeSize, binary, fileName );

    if ( ! common::Utils::validIndex( first, size ) )
    {
        array.clear();  
        return;        
    }

    IndexType nEntries = n;

    if ( n == nIndex )
    {
        nEntries = size - first;
    }
    else
    {
        SCAI_ASSERT_LE_ERROR( first + n, size, "array block size " << n << " invalid" )
    }

    // check if the specified data size fits the expected data type

    common::scalar::ScalarType dataType = mScalarTypeData;

    if ( mScalarTypeData == common::scalar::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    SCAI_ASSERT_EQUAL( dataTypeSize, (IndexType) common::typeSize( dataType ), 
                       "SAMG vector file has type size " << dataTypeSize 
                       << ", does not match to expected data type " << dataType )

    std::ios::openmode flags = std::ios::in;

    if ( binary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    IOStream inFile( dataFileName, flags );

    if ( binary )
    {
        inFile.skipBinary( first, dataTypeSize );
        inFile.readBinary( array, nEntries, dataType );
        inFile.skipBinary( size - nEntries - first, dataTypeSize );
    }
    else
    {
        inFile.skipFormatted( first );
        inFile.readFormatted( array, nEntries );
        inFile.skipFormatted( size - nEntries - first );
    }

    inFile.closeCheck();
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeMatrixHeader( 
    const IndexType numRows, 
    const IndexType numValues, 
    const bool binary, 
    const std::string& fileName )
{
    IOStream outFile( fileName, std::ios::out | std::ios::trunc );

    char fileType = binary ? 'b' : 'f';

    IndexType size = 1;
    IndexType rank = 0;

    outFile << fileType;
    outFile << " \t" << SAMG_IVERSION << "\n";
    outFile << "\t\t" << numValues;
    outFile << "\t" << numRows;
    outFile << "\t" << SAMG_VERSION_ID;
    outFile << "\t" << size;
    outFile << "\t" << rank;
    outFile.close();
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeStorageImpl(
    const MatrixStorage<ValueType>& storage,
    const std::string& fileName )
{
    utilskernel::LArray<IndexType> csrIA;
    utilskernel::LArray<IndexType> csrJA;
    utilskernel::LArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    // SAMG format starts indexing with 1

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    csrIA += 1;    
    csrJA += 1;     

    bool binary = ( mFileMode != FORMATTED );

    writeMatrixHeader( numRows, numValues, binary, fileName );

    SCAI_LOG_INFO( logger, *this << ": writeCSRData( " << fileName << " )" << ", #rows = " << csrIA.size() - 1
                            << ", #values = " << csrJA.size() )

    std::ios::openmode flags = std::ios::out | std::ios::trunc;

    if ( binary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    IOStream outFile( dataFileName, flags );

    if ( binary )
    {
        // take care of file type conversions as specified

        outFile.writeBinary( csrIA, mScalarTypeIndex );
        outFile.writeBinary( csrJA, mScalarTypeIndex ); 

        if ( mScalarTypeData != common::scalar::PATTERN )
        {
            outFile.writeBinary( csrValues, mScalarTypeData );
        }
    }
    else
    {
        int precIndex = 0;  
        int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

        // no conversions for formmatted write, but take care of precision

        outFile.writeFormatted( csrIA, precIndex );
        outFile.writeFormatted( csrJA, precIndex );

        if ( mScalarTypeData != common::scalar::PATTERN )
        {
            outFile.writeFormatted( csrValues, precData );
        }
    }

    outFile.close();
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readMatrixHeader( IndexType& numRows, IndexType& numValues, bool& binary, const std::string& fileName )
{
    IOStream inFile( fileName, std::ios::in );

    int  iversion; 
    char fileType = '!'; 

    IndexType id;
    IndexType size;
    IndexType rank;

    numValues = 0;

    inFile >> fileType >> iversion;

    if ( fileType == 'b' )
    {
        binary = true;
    }
    else if ( fileType == 'f' )
    {
        binary = false;
    }
    else 
    {   
        COMMON_THROWEXCEPTION( "Invalid file type = " << fileType << " in  SAMG header file " 
                                << fileName << ", must be either f or b" )
    }

    SCAI_ASSERT_EQUAL( SAMG_IVERSION, iversion, "SAMG version mismatch in SAMG file" << fileName )

    inFile >> numValues;
    inFile >> numRows;

    inFile >> id;
    inFile >> size;
    inFile >> rank; 

    inFile.close(); // explicitly, otherwise done by destructor

    SCAI_LOG_DEBUG( logger, "Info from header file " << fileName << ": #rows = " << numRows
                           << ", #values = " << numValues << ", type = " << fileType  )
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues, const std::string& fileName )
{
    bool binary;    // header only decides about formatted/binary read

    // start with reading the header

    readMatrixHeader( numRows, numValues, binary, fileName );

    numColumns = numRows;   // SAMG assumes always square matrices
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readStorageImpl(
    MatrixStorage<ValueType>& storage,
    const std::string& fileName,
    const IndexType firstRow,
    const IndexType nRows ) 
{
    IndexType numRows;   // SAMG always assumes square matrices
    IndexType numValues; // number of non-zero entries

    bool      binary;    // header only decides about formatted/binary read

    // start with reading the header

    readMatrixHeader( numRows, numValues, binary, fileName );

    SCAI_LOG_INFO( logger, "Info from header file " << fileName << ": #rows = " << numRows
                           << ", #values = " << numValues << ", binary = " << binary )

    if ( !common::Utils::validIndex( firstRow, numRows ) )
    {
        storage.clear();
        return;
    }

    IndexType numBlockRows = nRows;

    if ( nRows == nIndex )
    {
        numBlockRows = numRows - firstRow;
    }
    else
    {
        SCAI_ASSERT_LE_ERROR( firstRow + nRows, numRows, "storage block size " << numRows << " invalid" )
    }

    // now open the associated data file in correct mode

    std::ios::openmode flags = std::ios::in;

    if ( binary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    IOStream inFile( dataFileName, flags );

    utilskernel::LArray<IndexType> csrIA;
    utilskernel::LArray<IndexType> csrJA;
    utilskernel::LArray<ValueType> csrValues;

    size_t indexTypeSize = common::typeSize( mScalarTypeIndex );
    size_t valueTypeSize = sizeof( ValueType );

    if ( mScalarTypeData != common::scalar::INTERNAL )
    {
        valueTypeSize = common::typeSize( mScalarTypeData );
    }

    if ( binary )
    {
        // compare expected size with real size and give a warning

        size_t expectedSize = ( numRows + 1 + numValues ) * indexTypeSize + numValues * valueTypeSize;

        SCAI_LOG_INFO( logger, "expected size = " << expectedSize << ", type size = " << valueTypeSize 
                               << ", index size = " << indexTypeSize )

        inFile.seekg( 0, std::ios::end );
        size_t realSize = inFile.tellg();
        inFile.seekg( 0, std::ios::beg );

        if ( expectedSize != realSize )
        {
            SCAI_LOG_WARN( logger, "Binary file " << fileName << ": real size = " << realSize << 
                                    ", expected size = " << expectedSize << ", #rows = " << numRows << ", #nnz = " << numValues <<
                                    ", IndexType = " << mScalarTypeIndex << ", DataType = " << mScalarTypeData <<
                                    ", ValueType = " << common::TypeTraits<ValueType>::id() );
        }
    }

    if ( binary )
    {
        // Note: read operations can deal with scalar::INTERNAL, scalar::INDEX_TYPE

        inFile.skipBinary( firstRow, indexTypeSize );
        inFile.readBinary( csrIA, numBlockRows + 1, mScalarTypeIndex );
        inFile.skipBinary( numRows - numBlockRows - firstRow, indexTypeSize );
    }
    else
    {
        inFile.skipFormatted( firstRow );
        inFile.readFormatted( csrIA, numBlockRows + 1 );
        inFile.skipFormatted( numRows - numBlockRows - firstRow );
    }

    IndexType offsetBlock = csrIA[0];

    csrIA -= offsetBlock;    // offset array will now start at 0

    offsetBlock -= IndexType( 1 );        // SAMG indexing starts with 1, deal correctly for reading ja, values

    IndexType numBlockValues   = csrIA[numBlockRows];
 
    SCAI_LOG_DEBUG( logger, "ia        : read, offset = " << firstRow << ", numBlockRows = " << numBlockRows << " of " << numRows )
    SCAI_LOG_DEBUG( logger, "ja, values: read, offset = " << offsetBlock << ", numBlockValues = " << numBlockValues << " of " << numValues )

    if ( binary )
    {
        inFile.skipBinary( offsetBlock, indexTypeSize );
        inFile.readBinary( csrJA, numBlockValues, mScalarTypeIndex );
        inFile.skipBinary( numValues - offsetBlock - numBlockValues, indexTypeSize );
    }
    else
    {
        inFile.skipFormatted( offsetBlock );
        inFile.readFormatted( csrJA, numBlockValues );
        inFile.skipFormatted( numValues - offsetBlock - numBlockValues );
    }

    IndexType maxColumn = csrJA.max();   // maximal column index used

    csrJA -= IndexType( 1 );

    if ( mScalarTypeData == common::scalar::PATTERN )
    {
        csrValues.init( ValueType( 1 ), numBlockValues );   // set values with default value
    }
    else if ( binary )
    {
        inFile.skipBinary( offsetBlock, valueTypeSize );
        inFile.readBinary( csrValues, numBlockValues, mScalarTypeData );
        inFile.skipBinary( numValues - numBlockValues - offsetBlock, valueTypeSize );
    }
    else
    {
        inFile.skipFormatted( offsetBlock );
        inFile.readFormatted( csrValues, numBlockValues );
        inFile.skipFormatted( numValues - numBlockValues - offsetBlock );
    }

    inFile.closeCheck();   // gives a warning if not complete file has been read

    SCAI_LOG_INFO( logger, "CSR data: ia = " << csrIA << ", ja = " << csrJA << ", valaues = " << csrValues )

    IndexType numColumns = numRows;  // Usuallly, SAMG expects always square matrices
    
    if ( maxColumn > numColumns )
    {
        numColumns = maxColumn;      // but might be bigger for partitioned data
    }
 
    storage.setCSRData( numBlockRows, numColumns, numBlockValues, csrIA, csrJA, csrValues );
}

/* --------------------------------------------------------------------------------- */

SAMGIO::Guard SAMGIO::mGuard;

SAMGIO::Guard::Guard()
{
    addCreator( SAMG_VEC_HEADER_SUFFIX, &SAMGIO::create );
}

SAMGIO::Guard::~Guard()
{
    removeCreator( SAMG_VEC_HEADER_SUFFIX );
}

/* --------------------------------------------------------------------------------- */

int SAMGIO::deleteFile( const std::string& fileName )
{
    int rc = -1;

    if ( FileIO::hasSuffix( fileName, this->getMatrixFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else if ( FileIO::hasSuffix( fileName, this->getVectorFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else
    {
        SCAI_LOG_WARN( logger, *this << ", unsupported suffix for file " << fileName )
    }

    if ( rc != 0 ) 
    {
        return rc;
    }

    std::string dataFileName = getDataFileName( fileName );

    rc = std::remove( dataFileName.c_str() );

    return rc;
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
