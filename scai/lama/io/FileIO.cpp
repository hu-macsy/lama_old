/**
 * @file FileIO.cpp
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
 * @brief Implementation of non-pure methods for base class FileIO
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#include <scai/lama/io/FileIO.hpp>

#include <scai/lama/io/IOWrapper.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/exception/IOException.hpp>

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <string>
#include <fstream>
#include <memory>

namespace scai
{

namespace lama
{

std::ostream& operator<<( std::ostream& stream, const FileMode& value )
{
    switch ( value ) 
    {
        case FileMode::BINARY: 
            stream << "BINARY"; 
            break;
        case FileMode::FORMATTED: 
            stream << "FORMATTED"; 
            break;
        case FileMode::DEFAULT: 
            stream << "DEFAULT"; 
            break;
        default:
            stream << "UnknownFileMode";
    }
    return stream;
}

std::ostream& operator<<( std::ostream& stream, const DistributedIOMode& value )
{
    switch ( value ) 
    {
        case DistributedIOMode::MASTER: 
            stream << "MASTER"; 
            break;
        case DistributedIOMode::INDEPENDENT: 
            stream << "INDEPENDENT"; 
            break;
        case DistributedIOMode::COLLECTIVE: 
            stream << "COLLECTIVE"; 
            break;
        case DistributedIOMode::DEFAULT: 
            stream << "DEFAULT"; 
            break;
        default:
            stream << "UnknownDistributedIOMode";
    }
    return stream;
}

FileIO::FileIO() :

    mFileMode( FileMode::DEFAULT ),                      // no force of anything
    mDistMode( DistributedIOMode::DEFAULT ),             // no force of anything
    mScalarTypeIndex( common::ScalarType::INDEX_TYPE ),  // default is as used in LAMA
    mScalarTypeData( common::ScalarType::INTERNAL )      // default is same type as used in output data structure
{
    bool binary;

    if ( common::Settings::getEnvironment( binary, "SCAI_IO_BINARY" ) )
    {
        if ( binary )
        {
            mFileMode = FileMode::BINARY;
        }
        else
        {
            mFileMode = FileMode::FORMATTED;
        }

        SCAI_LOG_INFO( logger, "File mode set by SCAI_IO_BINARY = " << binary << " : " << mFileMode )
    }

    std::string datatype;

    if ( common::Settings::getEnvironment( datatype, "SCAI_IO_TYPE_DATA" ) )
    {
        mScalarTypeData = common::str2ScalarType( datatype.c_str() );

        if ( common::ScalarType::UNKNOWN == mScalarTypeData )
        {
            COMMON_THROWEXCEPTION( "Not a known value type: SCAI_IO_TYPE_DATA="  << datatype )
        }
    }

    if ( common::Settings::getEnvironment( datatype, "SCAI_IO_TYPE_INDEX" ) )
    {
        mScalarTypeIndex = common::str2ScalarType( datatype.c_str() );

        if ( common::ScalarType::UNKNOWN == mScalarTypeIndex )
        {
            COMMON_THROWEXCEPTION( "Not a known value type: SCAI_IO_TYPE_INDEX="  << datatype )
        }
    }

    SCAI_LOG_DEBUG( logger, "Constructed FileIO" )
}

/* --------------------------------------------------------------------------------- */

FileIO::~FileIO()
{
}

/* --------------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( FileIO::logger, "FileIO" )

/* --------------------------------------------------------------------------------- */

bool FileIO::hasCollectiveIO() const
{
    // collective I/O is not supported by default

    return false;
}

/* --------------------------------------------------------------------------------- */

void FileIO::open( const char* fileName, const char* openMode, const DistributedIOMode distMode )
{
    SCAI_LOG_INFO( logger, *this << ": open, name = " << fileName << ", open mode = " << openMode )

    mDistMode = distMode;

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    if ( strstr( fileName, "%r" ) != NULL )
    {
        mDistMode = DistributedIOMode::INDEPENDENT;

        std::string independentFileName = fileName;

        size_t pos = independentFileName.find( "%r" );

        std::ostringstream rankStr;

        PartitionId size = comm->getSize();
        PartitionId rank = comm->getRank();

        if ( size > 1 )
        {
            rankStr << rank << "." << size;
        }

        independentFileName.replace( pos, 2, rankStr.str() );

        openIt( independentFileName.c_str(), openMode );
    }
    else
    {
        if ( mDistMode == DistributedIOMode::DEFAULT )
        {
            if ( hasCollectiveIO() )
            {
                mDistMode = DistributedIOMode::COLLECTIVE;
            }
            else
            {
                mDistMode = DistributedIOMode::MASTER;
            }
        }

        // in single mode only MASTER processor opens it

        if ( mDistMode != DistributedIOMode::MASTER || comm->getRank() == 0 )
        {
            openIt( fileName, openMode );
        }
    }

    SCAI_LOG_DEBUG( logger, *this << " is now open" )
}

/* --------------------------------------------------------------------------------- */

void FileIO::close()
{
    auto comm = getCommunicatorPtr();

    // in single mode only MASTER processor closes it

    if ( mDistMode != DistributedIOMode::MASTER || comm->getRank() == 0 )
    {
        closeIt();
    }

    // synchronize at close avoids that file is seen by all processors in same state

    if ( mDistMode == DistributedIOMode::MASTER )
    {
        comm->synchronize();
    }

    mDistMode = DistributedIOMode::DEFAULT;
}

/* --------------------------------------------------------------------------------- */

void FileIO::writeAt( std::ostream& stream ) const
{
    stream << "FileIO( ";
    writeMode( stream );
    stream << " )";
}

/* --------------------------------------------------------------------------------- */

void FileIO::writeMode( std::ostream& stream ) const
{
    stream << "FileMode = " << mFileMode;
    stream << ", DistributedIOMode = " << mDistMode;
    stream << ", data = " << mScalarTypeData;
    stream << ", index = " << mScalarTypeIndex;
}

/* --------------------------------------------------------------------------------- */

int FileIO::deleteFile( const std::string& fileName )
{
    int rc = -1;

    if ( hasSuffix( fileName, getMatrixFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else if ( hasSuffix( fileName, getVectorFileSuffix() ) )
    {
        rc = std::remove( fileName.c_str() );
    }
    else
    {
        SCAI_LOG_WARN( logger, "do not delete file with unknown suffix" )
    }

    return rc;
}

/* --------------------------------------------------------------------------------- */

void FileIO::setIndexType( common::ScalarType type )
{
    mScalarTypeIndex = type;
}

void FileIO::setDataType( common::ScalarType type )
{
    mScalarTypeData = type;
}

void FileIO::setMode( const FileMode mode )
{
    // only a forced file mode overwrites the current file mode

    if ( mode != FileMode::DEFAULT )
    {
        mFileMode = mode;
    }
}

/* --------------------------------------------------------------------------------- */

int FileIO::getDataPrecision( common::ScalarType stype )
{
    int prec = 0;

    if ( mScalarTypeData == common::ScalarType::INTERNAL )
    {
        prec = common::precision( stype );
    }
    else
    {
        prec = common::precision( mScalarTypeData );
    }

    // user can have redefined the precision

    common::Settings::getEnvironment( prec, "SCAI_IO_PRECISION" );

    return prec;
}

/* --------------------------------------------------------------------------------- */

bool FileIO::fileExists( const std::string& fileName )
{
    // open the file for reading -> it exists
    std::ifstream file( fileName.c_str(), std::ios::in );
    return file.good();
}

/* -------------------------------------------------------------------------- */

bool FileIO::hasSuffix( const std::string& fileName, const std::string& suffix )
{
    size_t suffixSize = suffix.size();

    // Note: hasSuffix( ".frv", ".frv" ) returns also false, avoids empty names

    return fileName.size() > suffixSize &&
           fileName.compare( fileName.size() - suffixSize, suffixSize, suffix ) == 0;
}

/* -------------------------------------------------------------------------- */

std::string FileIO::getSuffix( const std::string& fileName )
{
    size_t pos = fileName.find_last_of( "." );

    if ( pos == std::string::npos )
    {
        return "";
    }

    return fileName.substr( pos );
}

/* -------------------------------------------------------------------------- */

int FileIO::removeFile( const std::string& fileName )
{
    std::string suffix = getSuffix( fileName );

    // Let FileIO class delete the file as there might be joint files

    if ( canCreate( suffix ) )
    {
        std::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        SCAI_LOG_INFO( logger, "delete file via " << *fileIO )

        return fileIO->deleteFile( fileName );
    }

    // unknwon suffix, also delete it

    int rc = std::remove( fileName.c_str() );

    return rc;
}

/* -------------------------------------------------------------------------- */

void FileIO::write(
    const hmemo::_HArray& array,
    const std::string& outFileName,
    const common::ScalarType dataType )
{
    // This static method creates a FileIO object by suffix (file type)  and calls its write method for dense array

    std::string suffix = getSuffix( outFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "Unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    fileIO->setDataType( dataType );
    fileIO->open( outFileName.c_str(), "w", DistributedIOMode::INDEPENDENT );
    fileIO->writeArray( array );
    fileIO->close();
}

/* -------------------------------------------------------------------------- */

void FileIO::write(
    const IndexType size,
    const void* zero,
    const hmemo::HArray<IndexType>& indexes,
    const hmemo::_HArray& values,
    const std::string& outFileName,
    const common::ScalarType dataType )
{
    // This static method creates a FileIO object by suffix (file type)  and calls its write method for sparse array

    std::string suffix = getSuffix( outFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "Unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    fileIO->setDataType( dataType );
    fileIO->open( outFileName.c_str(), "w", DistributedIOMode::INDEPENDENT );
    fileIO->writeSparse( size, zero, indexes, values );
    fileIO->close();
}

/* -------------------------------------------------------------------------- */

void FileIO::read(
    hmemo::_HArray& array,
    const std::string& inFileName,
    const common::ScalarType dataType )
{
    std::string suffix = getSuffix( inFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "ERROR: read from file " << inFileName <<
                             ": unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    fileIO->setDataType( dataType );
    fileIO->open( inFileName.c_str(), "r", DistributedIOMode::INDEPENDENT );
    fileIO->readArray( array );
    fileIO->close();
}

/* -------------------------------------------------------------------------- */

void FileIO::read(
    hmemo::_HArray& array,
    common::Grid& grid,
    const std::string& inFileName,
    const common::ScalarType dataType )
{
    std::string suffix = getSuffix( inFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "ERROR: read from file " << inFileName <<
                             ": unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    fileIO->setDataType( dataType );
    fileIO->open( inFileName.c_str(), "r", DistributedIOMode::INDEPENDENT );
    fileIO->readGridArray( array, grid );
    fileIO->close();
}

/* -------------------------------------------------------------------------- */

void FileIO::read(
    IndexType& size,
    void* zero,
    hmemo::HArray<IndexType>& indexes,
    hmemo::_HArray& values,
    const std::string& inFileName,
    const common::ScalarType dataType )
{
    std::string suffix = getSuffix( inFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "ERROR: read from file " << inFileName <<
                             ": unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    fileIO->setDataType( dataType );
    fileIO->open( inFileName.c_str(), "r", DistributedIOMode::INDEPENDENT );
    fileIO->readSparse( size, zero, indexes, values );
    fileIO->close();
}

/* -------------------------------------------------------------------------- */

IndexType FileIO::getArraySize( const std::string& inFileName )
{
    std::string suffix = getSuffix( inFileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "ERROR: read from file " << inFileName <<
                             ": unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    IndexType size = invalidIndex;

    fileIO->open( inFileName.c_str(), "r", DistributedIOMode::INDEPENDENT );
    fileIO->getArrayInfo( size );
    fileIO->close();

    return size;
}

/* -------------------------------------------------------------------------- */

IndexType FileIO::getStorageSize( const std::string& fileName )
{
    std::string suffix = getSuffix( fileName );

    if ( !canCreate( suffix ) )
    {
        SCAI_THROWEXCEPTION( common::IOException, "Unsupported suffix " << suffix << ", no FileIO handler availabe" )
    }

    std::unique_ptr<FileIO> fileIO ( FileIO::create( suffix ) );

    IndexType size = invalidIndex;

    IndexType numColumns = 0;   // dummy
    IndexType numValues = 0;    // dummy

    fileIO->open( fileName.c_str(), "r", DistributedIOMode::INDEPENDENT );
    fileIO->getStorageInfo( size, numColumns, numValues );
    fileIO->close();

    return size;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void FileIO::writeSparseImpl(
    const IndexType size,
    const ValueType& zero,
    const hmemo::HArray<IndexType>& indexes,
    const hmemo::HArray<ValueType>& values )
{
    // sparse unsupported for this file format, write it dense

    hmemo::HArray<ValueType> denseArray;
    utilskernel::HArrayUtils::buildDenseArray( denseArray, size, values, indexes, zero );
    writeArray( denseArray );
}

void FileIO::writeSparse( const IndexType n, const void* zero, const hmemo::HArray<IndexType>& indexes, const hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<FileIO, SCAI_ARRAY_TYPES_HOST_LIST>::writeSparse( *this, n, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void FileIO::readSparseImpl(
    IndexType& size,
    ValueType& zero,
    hmemo::HArray<IndexType>& indexes,
    hmemo::HArray<ValueType>& values )
{
    // sparse array not supported for this file format, read a dense array

    SCAI_LOG_WARN( logger, "read sparse array not supported for " << *this << ", reads dense array." )

    hmemo::HArray<ValueType> denseArray;

    readArray( denseArray );
    size = denseArray.size();
    zero = 0;
    utilskernel::HArrayUtils::buildSparseArray( values, indexes, denseArray, zero );
}

void FileIO::readSparse( IndexType& size, void* zero, hmemo::HArray<IndexType>& indexes, hmemo::_HArray& values )
{
    // use IOWrapper to called the typed version of this routine

    IOWrapper<FileIO, SCAI_ARRAY_TYPES_HOST_LIST>::readSparse( *this, size, zero, indexes, values );
}

/* --------------------------------------------------------------------------------- */

DistributedIOMode FileIO::getDistributedIOMode() const
{
    // mDistMode == DEFAULT indicates that no file is open at all

    return mDistMode;
}

/* --------------------------------------------------------------------------------- */

dmemo::CommunicatorPtr FileIO::getCommunicatorPtr() const
{
    return dmemo::Communicator::getCommunicatorPtr();
}

/* --------------------------------------------------------------------------------- */


}  // namespace lama

}  // namespace scai
