/**
 * @file FileIO.cpp
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
 * @brief Implementation of non-pure methods for base class FileIO
 * @author Thomas Brandes
 * @date 20.06.2016
 */

#include "FileIO.hpp"

#include <scai/common/Settings.hpp>
#include <string>
#include <fstream>

namespace scai
{

namespace lama
{

FileIO::FileIO() :

    mFileMode( DEFAULT_MODE ),                       // no forace of anything
    mAppendMode( false ),                            // default is to write each output file new
    mScalarTypeIndex( common::scalar::INDEX_TYPE ),  // default is as used in LAMA
    mScalarTypeData( common::scalar::INTERNAL )      // default is same type as used in output data structure
{
    bool binary;

    if ( common::Settings::getEnvironment( binary, "SCAI_IO_BINARY" ) )
    {
        if ( binary )
        {
            mFileMode = BINARY;
        }
        else
        {
            mFileMode = FORMATTED;
        }

        SCAI_LOG_INFO( logger, "File mode set by SCAI_IO_BINARY = " << binary )
    }

    common::Settings::getEnvironment( mAppendMode, "SCAI_IO_APPEND" );

    std::string datatype;

    if ( common::Settings::getEnvironment( datatype, "SCAI_IO_TYPE_DATA" ) )
    {
        mScalarTypeData = common::str2ScalarType( datatype.c_str() );
        
        if ( common::scalar::UNKNOWN == mScalarTypeData )
        {
            COMMON_THROWEXCEPTION( "Not a known value type: SCAI_IO_TYPE_DATA="  << datatype )
        }
    }

    if ( common::Settings::getEnvironment( datatype, "SCAI_IO_TYPE_INDEX" ) )
    {
        mScalarTypeIndex = common::str2ScalarType( datatype.c_str() );

        if ( common::scalar::UNKNOWN == mScalarTypeIndex )
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

void FileIO::writeAt( std::ostream& stream ) const
{
    stream << "FileIO( ";
    writeMode( stream );
    stream << " )";
}

/* --------------------------------------------------------------------------------- */

void FileIO::writeMode( std::ostream& stream ) const
{
    stream << "FileMode = ";

    if ( mFileMode == BINARY )
    {
        stream << "binary";
    }
    else if ( mFileMode == FORMATTED )
    {
        stream << "formatted";
    }
    else
    {
        stream << "DEFAULT";
    }

    stream << ", append = " << mAppendMode;
    stream << ", data = " << mScalarTypeData;
    stream << ", index = " << mScalarTypeIndex;
}

/* --------------------------------------------------------------------------------- */

void FileIO::setIndexType( common::scalar::ScalarType type ) 
{
    mScalarTypeIndex = type;
}

void FileIO::setDataType( common::scalar::ScalarType type ) 
{
    mScalarTypeData = type;
}

void FileIO::setMode( const FileMode mode )
{
    mFileMode = mode;
}

void FileIO::enableAppendMode( bool flag ) 
{
    mAppendMode = flag;
}

/* --------------------------------------------------------------------------------- */

int FileIO::getDataPrecision( common::scalar::ScalarType stype )
{
    int prec = 0;

    if ( mScalarTypeData == common::scalar::INTERNAL )
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
        common::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        return fileIO->deleteFile( fileName );
    }

    // unknwon suffix, also delete it

    int rc = std::remove( fileName.c_str() );

    return rc;
}

}  // namespace lama

}  // namespace scai
