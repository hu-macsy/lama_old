/**
 * @file SAMGIO.cpp
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
 * @brief Implementation of IO methods for SAMG format
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#include <scai/lama/io/SAMGIO.hpp>

#include <scai/lama/io/IOStream.hpp>
#include <scai/lama/io/IOWrapper.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

#include <scai/common/Grid.hpp>
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
using utilskernel::HArrayUtils;

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

    if ( FileIO::hasSuffix( headerFileName, SAMG_MAT_HEADER_SUFFIX ) )
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

    SCAI_LOG_INFO( logger, "SAMGIO default object: " << *this )
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::openIt( const std::string& fileName, const char* openMode )
{
    SCAI_LOG_INFO( logger, "SAMGIO: openIt ( fileName = " << fileName << ", openMode = " << openMode << " ), fileMode = " << mFileMode )
                           
    std::ios::openmode flags;

    if ( strcmp( openMode, "w" ) == 0 )
    {
        flags = std::ios::out | std::ios::trunc;
    }
    else if ( strcmp( openMode, "r" ) == 0 )
    {
        flags = std::ios::in;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported file mode for SAMG file: " << openMode )
    }

    mHeaderFile.open( fileName, flags );

    if ( strcmp( openMode, "w" ) == 0 )
    {
        mBinary = mFileMode != FileMode::FORMATTED;
    }
    else if ( strcmp( openMode, "r" ) == 0 )
    {
        readFileMode();
        mHeaderFile.seekg( 0, std::ios::beg );
    }

    if ( mBinary )
    {
        flags |= std::ios::binary;
    }

    std::string dataFileName = getDataFileName( fileName );

    mDataFile.open( dataFileName, flags );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readFileMode()
{
    char fileType = ' ';

    mHeaderFile >> fileType;

    if ( fileType == 'b' )
    {
        mBinary = true;
    }
    else if ( fileType == 'f' )
    {
        mBinary = false;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << mHeaderFile.getFileName()
                               << ", fileType = " << fileType << " illegal" )
    }

    // ToDo: verify if file mode matches dataFile if it has already been opened
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::closeIt()
{
    mHeaderFile.close();
    mDataFile.close();
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeVectorHeader(
    const IndexType n,
    const IndexType typeSize )
{
    char fileType = mBinary ? 'b' : 'f';

    mHeaderFile << fileType << std::endl;
    mHeaderFile << n << std::endl;
    mHeaderFile << typeSize;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeArrayImpl( const HArray<ValueType>& array )
{

    SCAI_LOG_INFO( logger, "writeArrayImpl<" << common::TypeTraits<ValueType>::id() << ">, array = " << array
                   << " to " << mHeaderFile.getFileName() << ", " << mDataFile.getFileName() )

    // needed for header file: type size is size of data type used in output

    int typeSize = sizeof( ValueType );

    if ( mScalarTypeData != common::ScalarType::INTERNAL )
    {
        typeSize = common::typeSize( mScalarTypeData );
    }

    writeVectorHeader( array.size(), typeSize );

    if ( mBinary )
    {
        mDataFile.writeBinary( array, mScalarTypeData );
    }
    else
    {
        int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

        mDataFile.writeFormatted( array, precData );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeSparseImpl(
    const IndexType size,
    const ValueType& zero,
    const HArray<IndexType>& indexes,
    const HArray<ValueType>& values )
{
    // sparse unsupported for SAMG file format, write it dense

    SCAI_LOG_INFO( logger, "writeSparseImpl, size = " << size << ", nnz = " << values.size()
                   << " to " << mHeaderFile.getFileName() );

    HArray<ValueType> denseArray;
    utilskernel::HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );

    writeArrayImpl( denseArray );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readVectorHeader( IndexType& n, IndexType& typeSize )
{
    typeSize = 0;
    n        = 0;

    readFileMode();

    mHeaderFile >> n;
    mHeaderFile >> typeSize;

    if ( mHeaderFile.fail() )
    {
        COMMON_THROWEXCEPTION( "Invalid SAMG vector header file: " << mHeaderFile.getFileName()
                               << ", could not read '[f|b] <n> <typeSize>" )
    }
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::getArrayInfo( IndexType& size )
{
    IndexType dataTypeSize;   // dummy variable needed for readVectorHeader

    std::streampos pos = mHeaderFile.tellg();

    readVectorHeader( size, dataTypeSize );

    SCAI_LOG_INFO( logger, "read  array info, size = " << size << ", data type size = " << dataTypeSize )

    mHeaderFile.clear();       // important to reset flags
    mHeaderFile.seekg( pos );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readArray( hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<SAMGIO, SCAI_ARRAY_TYPES_HOST_LIST>::readArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readArrayImpl( HArray<ValueType>& array )
{
    IndexType dataTypeSize = 0;
    IndexType size = 0;

    // start with reading the *.frv header file

    readVectorHeader( size, dataTypeSize );

    IndexType nEntries = size;

    // check if the specified data size fits the expected data type

    common::ScalarType dataType = mScalarTypeData;

    if ( mScalarTypeData == common::ScalarType::INTERNAL )
    {
        dataType = common::TypeTraits<ValueType>::stype;
    }

    SCAI_ASSERT_EQUAL( dataTypeSize, ( IndexType ) common::typeSize( dataType ),
                       "SAMG vector file has type size " << dataTypeSize
                       << ", does not match to expected data type " << dataType )

    if ( mBinary )
    {
        mDataFile.readBinary( array, nEntries, dataType );
    }
    else
    {
        mDataFile.readFormatted( array, nEntries );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readSparseImpl(
    IndexType& size,
    ValueType& zero,
    HArray<IndexType>& indexes,
    HArray<ValueType>& values )
{
    // sparse array not supported for this file format, uses a temporary dense array of same type

    HArray<ValueType> denseArray;

    readArray( denseArray );
    size = denseArray.size();
    zero = 0;
    utilskernel::HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeMatrixHeader(
    const IndexType numRows,
    const IndexType numValues )
{
    char fileType = mBinary ? 'b' : 'f';

    IndexType size = 1;
    IndexType rank = 0;

    mHeaderFile << fileType;
    mHeaderFile << " \t" << SAMG_IVERSION << "\n";
    mHeaderFile << "\t\t" << numValues;
    mHeaderFile << "\t" << numRows;
    mHeaderFile << "\t" << SAMG_VERSION_ID;
    mHeaderFile << "\t" << size;
    mHeaderFile << "\t" << rank;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::writeStorageImpl( const MatrixStorage<ValueType>& storage )
{
    SCAI_LOG_INFO( logger, "write storage ( mBinary = " << mBinary << " ): " << storage )

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    storage.buildCSRData( csrIA, csrJA, csrValues );

    // SAMG format starts indexing with 1

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    HArrayUtils::compute<IndexType>( csrIA, csrIA, common::BinaryOp::ADD, 1 );
    HArrayUtils::compute<IndexType>( csrJA, csrJA, common::BinaryOp::ADD, 1 );

    writeMatrixHeader( numRows, numValues );

    SCAI_LOG_INFO( logger, *this << ": writeCSRData( " << mHeaderFile.getFileName() << " )" << ", #rows = " << csrIA.size() - 1
                   << ", #values = " << csrJA.size() )

    if ( mBinary )
    {
        // take care of file type conversions as specified

        mDataFile.writeBinary( csrIA, mScalarTypeIndex );
        mDataFile.writeBinary( csrJA, mScalarTypeIndex );

        if ( mScalarTypeData != common::ScalarType::PATTERN )
        {
            mDataFile.writeBinary( csrValues, mScalarTypeData );
        }
    }
    else
    {
        int precIndex = 0;
        int precData  = getDataPrecision( common::TypeTraits<ValueType>::stype );

        // no conversions for formmatted write, but take care of precision

        mDataFile.writeFormatted( csrIA, precIndex );
        mDataFile.writeFormatted( csrJA, precIndex );

        if ( mScalarTypeData != common::ScalarType::PATTERN )
        {
            mDataFile.writeFormatted( csrValues, precData );
        }
    }
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readMatrixHeader( IndexType& numRows, IndexType& numValues )
{
    int  iversion;
    char fileType = '!';

    IndexType id;
    IndexType size;
    IndexType rank;

    numValues = 0;

    readFileMode();

    const std::string& fileName = mHeaderFile.getFileName();

    mHeaderFile >> iversion;

    SCAI_ASSERT_EQUAL( SAMG_IVERSION, iversion, "SAMG version mismatch in SAMG file" << fileName )

    mHeaderFile >> numValues;
    mHeaderFile >> numRows;

    mHeaderFile >> id;
    mHeaderFile >> size;
    mHeaderFile >> rank;

    SCAI_LOG_DEBUG( logger, "Info from header file " << fileName << ": #rows = " << numRows
                    << ", #values = " << numValues << ", type = " << fileType  )
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::getStorageInfo( IndexType& numRows, IndexType& numColumns, IndexType& numValues )
{
    // start with reading the header

    std::streampos pos = mHeaderFile.tellg();

    readMatrixHeader( numRows, numValues );

    mHeaderFile.clear();       // important to reset flags
    mHeaderFile.seekg( pos );

    numColumns = numRows;   // SAMG assumes always square matrices
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void SAMGIO::readStorageImpl( MatrixStorage<ValueType>& storage )
{
    IndexType numRows;   // SAMG always assumes square matrices
    IndexType numValues; // number of non-zero entries

    // start with reading the header

    readMatrixHeader( numRows, numValues );

    SCAI_LOG_INFO( logger, "Info from header file " << mHeaderFile.getFileName() << ": #rows = " << numRows
                   << ", #values = " << numValues << ", binary = " << mBinary )

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    size_t indexTypeSize = common::typeSize( mScalarTypeIndex );
    size_t valueTypeSize = sizeof( ValueType );

    if ( mScalarTypeData != common::ScalarType::INTERNAL )
    {
        valueTypeSize = common::typeSize( mScalarTypeData );
    }

    if ( mBinary )
    {
        // compare expected size with real size and give a warning

        size_t expectedSize = ( numRows + 1 + numValues ) * indexTypeSize + numValues * valueTypeSize;

        SCAI_LOG_INFO( logger, "expected size = " << expectedSize << ", type size = " << valueTypeSize
                       << ", index size = " << indexTypeSize )

        mDataFile.seekg( 0, std::ios::end );
        size_t realSize = mDataFile.tellg();
        mDataFile.seekg( 0, std::ios::beg );

        if ( expectedSize != realSize )
        {
            SCAI_LOG_WARN( logger, "Binary file " << mDataFile.getFileName() << ": real size = " << realSize <<
                           ", expected size = " << expectedSize << ", #rows = " << numRows << ", #nnz = " << numValues <<
                           ", IndexType = " << mScalarTypeIndex << ", DataType = " << mScalarTypeData <<
                           ", ValueType = " << common::TypeTraits<ValueType>::id() );
        }
    }

    if ( mBinary )
    {
        // Note: read operations can deal with ScalarType::INTERNAL, ScalarType::INDEX_TYPE

        mDataFile.readBinary( csrIA, numRows + 1, mScalarTypeIndex );
    }
    else
    {
        mDataFile.readFormatted( csrIA, numRows + 1 );
    }

    IndexType offset = csrIA[0];

    HArrayUtils::compute<IndexType>( csrIA, csrIA, common::BinaryOp::SUB, offset );   // offset array will now start at 0

    if ( mBinary )
    {
        mDataFile.readBinary( csrJA, numValues, mScalarTypeIndex );
    }
    else
    {
        mDataFile.readFormatted( csrJA, numValues );
    }

    IndexType maxColumn = utilskernel::HArrayUtils::max( csrJA );   // maximal column index used

    HArrayUtils::compute<IndexType>( csrJA, csrJA, common::BinaryOp::SUB, 1 );

    if ( mScalarTypeData == common::ScalarType::PATTERN )
    {
        csrValues.setSameValue( numValues, ValueType( 1 ) );   // set values with default value
    }
    else if ( mBinary )
    {
        mDataFile.readBinary( csrValues, numValues, mScalarTypeData );
    }
    else
    {
        mDataFile.readFormatted( csrValues, numValues );
    }

    mDataFile.closeCheck();   // gives a warning if not complete file has been read

    SCAI_LOG_INFO( logger, "CSR data: ia = " << csrIA << ", ja = " << csrJA << ", valaues = " << csrValues )

    IndexType numColumns = numRows;  // Usuallly, SAMG expects always square matrices

    if ( maxColumn > numColumns )
    {
        numColumns = maxColumn;      // but might be bigger for partitioned data
    }

    storage.setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeGridArray( const hmemo::_HArray& data, const common::Grid& grid )
{
    if ( grid.nDims() > 1 )
    {
        SCAI_LOG_WARN( logger, "Grid shape information is lost for array when writing to file" )
    }

    writeArray( data );
}

void SAMGIO::readGridArray( hmemo::_HArray& data, common::Grid& grid )
{
    readArray( data );
    grid = common::Grid1D( data.size() );
    SCAI_LOG_WARN( logger, "SAMG does not support multidimensional array, take default shape " << grid )
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

void SAMGIO::writeStorage( const _MatrixStorage& storage )
{
    IOWrapper<SAMGIO, SCAI_NUMERIC_TYPES_HOST_LIST>::writeStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readStorage( _MatrixStorage& storage )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<SAMGIO, SCAI_NUMERIC_TYPES_HOST_LIST>::readStorage( *this, storage );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeArray( const hmemo::_HArray& array )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<SAMGIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeArray( *this, array );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::writeSparse( const IndexType n, const void* zero, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<SAMGIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparse( *this, n, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

void SAMGIO::readSparse( IndexType& size, void* zero, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<SAMGIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparse( *this, size, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

}  // lama

}  // scai
