/**
 * @file NoCollectiveFile.cpp
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
 * @brief Collectice/concurrent file I/O using MPI I/O
 * @author Thomas Brandes
 * @date 14.01.2019
 */

#include <scai/dmemo/NoCollectiveFile.hpp>

#include <scai/common/exception/IOException.hpp>

namespace scai
{

namespace dmemo
{

/* ---------------------------------------------------------------------------------- */

NoCollectiveFile::NoCollectiveFile( CommunicatorPtr comm ) :

   CollectiveFile( comm )
{
}

/* ---------------------------------------------------------------------------------- */

NoCollectiveFile::~NoCollectiveFile()
{
    if ( mFile != NULL )
    {
        close();
    }
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::open( const char* fileName, const char* fileMode )
{
    mFile = fopen( fileName, fileMode );

    if ( mFile == NULL )
    {
        SCAI_THROWEXCEPTION( common::IOException,
                             "could not open file " << fileName << ", mode = " << fileMode )
    }
    else
    {
        CollectiveFile::set( fileName, 0 );
    }
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::close()
{
    if ( mFile == NULL )
    {
        SCAI_LOG_WARN( logger, "call of close, file already closed" )
    }
    else
    {
        fclose( mFile );
        mFile = NULL;
    }

    CollectiveFile::set( "", 0 );
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    SCAI_LOG_INFO( logger, "write " << n << " entries " << stype << " @ offset = " << offset )
    SCAI_ASSERT_ERROR( mFile, "no file open" )
    fseek( mFile, offset, SEEK_SET );
    // n == 0 : val might be null that would result in sanitize error
    if ( n > 0 ) 
    {
        size_t writeN = fwrite( val, typeSize( stype ), n, mFile );
        SCAI_ASSERT_EQ_ERROR( writeN, n, "insufficent write" )
    }
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    writeSingleImpl( offset, val, n, stype );
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    SCAI_LOG_INFO( logger, "read " << n << " entries " << stype << " @ offset = " << offset )

    SCAI_ASSERT_ERROR( mFile, "no file open" )
    fseek( mFile, offset, SEEK_SET );
    size_t elemSize = typeSize( stype );
    size_t readN = fread( val, elemSize, n, mFile );
    SCAI_ASSERT_EQ_ERROR( readN, n, "insufficent read" )
}

/* ---------------------------------------------------------------------------------- */

void NoCollectiveFile::readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    readSingleImpl( val, n, offset, stype );
}

}

}
