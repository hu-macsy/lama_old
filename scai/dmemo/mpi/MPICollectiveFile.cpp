/**
 * @file MPICollectiveFile.cpp
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

#include <scai/dmemo/mpi/MPICollectiveFile.hpp>

#include <scai/dmemo/mpi/MPIUtils.hpp>
#include <scai/tracing.hpp>

#include <scai/common/exception/IOException.hpp>

namespace scai
{

namespace dmemo
{

/* -------------------------------------------------------------------------- */

MPICollectiveFile::MPICollectiveFile( CommunicatorPtr comm ) :

   CollectiveFile( comm )

{
}

/* -------------------------------------------------------------------------- */

void MPICollectiveFile::open( const char* fileName, const char* fileMode )
{
    SCAI_REGION( "ColFile.MPI.open" )

    SCAI_LOG_INFO( logger, "MPI open collecitve file " << fileName << ", mode " << fileMode )

    // we need the MPI communicator from the communicator

    auto mpiComm = std::dynamic_pointer_cast<const MPICommunicator>( mComm );

    SCAI_ASSERT( mpiComm, "MPICollective file only with MPI communicator, comm = " << *mComm )

    int openMode = 0;

    if ( strcmp( fileMode, "r" ) == 0 )
    {
        openMode = MPI_MODE_RDONLY;
    }
    else if ( strcmp( fileMode, "w" ) == 0 )
    {
        openMode = MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_WRONLY;

        int rc = MPI_File_open( mpiComm->getMPIComm(), const_cast<char*>( fileName ), openMode, MPI_INFO_NULL, &mFileHandle );

        SCAI_LOG_DEBUG( logger, "open file " << fileName << " for delete, rc = " << rc )

        if ( rc == MPI_SUCCESS )
        {
            SCAI_MPICALL( logger, MPI_File_close( &mFileHandle ), "could not delete existing file " << fileName )
        }

        openMode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL;
    }
    else if ( strcmp( fileMode, "a" ) == 0 )
    {
        openMode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND;
    }
    else
    {
        SCAI_THROWEXCEPTION( common::IOException, "illegal mode = " << fileMode << " to open file " << fileName )
    }

    SCAI_MPICALL( logger,
                  MPI_File_open( mpiComm->getMPIComm(), const_cast<char*>( fileName ), openMode, MPI_INFO_NULL, &mFileHandle ),
                  "MPI_File_open: fileName = " << fileName << ", mode = " << fileMode )

    if ( strcmp( fileMode, "a" ) == 0 )
    {
        CollectiveFile::set( fileName, getSize() );
    }
    else
    {
        CollectiveFile::set( fileName, 0 );
    }
}

/* -------------------------------------------------------------------------- */

void MPICollectiveFile::close()
{
    SCAI_REGION( "ColFile.MPI.close" )

    SCAI_LOG_INFO( logger, "close file " << mFileName << ", offset = " << mOffset );

    SCAI_MPICALL( logger, MPI_File_close( &mFileHandle ), "MPI_File_close" )
}

/* -------------------------------------------------------------------------- */

size_t MPICollectiveFile::writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    SCAI_LOG_INFO( logger, *mComm << ": writeSingle( offset = " << offset << ", size = " << n << ", type = " << stype )

    MPI_Status stat;

    MPI_Datatype commType = MPICommunicator::getMPIType( stype );

    SCAI_MPICALL( logger, 
      MPI_File_write_at( mFileHandle, offset, const_cast<void*>( val ), n, commType, &stat ),
      "MPI_File_write_at" )

    int writtenEntries = 0;

    SCAI_MPICALL( logger, MPI_Get_count( &stat, commType, &writtenEntries ),
                  "MPI_Get_count<" << stype << ">" )

    return static_cast<size_t>( writtenEntries );
}

/* -------------------------------------------------------------------------- */

size_t MPICollectiveFile::writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    SCAI_LOG_INFO( logger, *mComm << ": writeAll( offset = " << offset << ", size = " << n << ", type = " << stype )

    MPI_Status stat;

    MPI_Datatype commType = MPICommunicator::getMPIType( stype );

    SCAI_MPICALL( logger,
        MPI_File_write_at_all(
            mFileHandle, 
            offset, 
            const_cast<void*>( val ), n,
            commType, &stat ),
        "MPI_File_write_at_all" )

    int writtenEntries = 0;

    SCAI_MPICALL( logger, MPI_Get_count( &stat, commType, &writtenEntries ),
                  "MPI_Get_count<" << stype << ">" )

    return static_cast<size_t>( writtenEntries );
}

/* -------------------------------------------------------------------------- */

size_t MPICollectiveFile::readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    MPI_Status stat;

    MPI_Datatype commType = MPICommunicator::getMPIType( stype );

    SCAI_MPICALL( logger,
        MPI_File_read_at( mFileHandle, offset, val, n, commType, &stat ),
        "MPI_File_read_at" )

    int readEntries = 0;

    SCAI_MPICALL( logger, MPI_Get_count( &stat, commType, &readEntries ),
                  "MPI_Get_count<" << stype << ">" )

    return static_cast<size_t>( readEntries );
}

/* -------------------------------------------------------------------------- */

size_t MPICollectiveFile::readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    MPI_Status stat;

    MPI_Datatype commType = MPICommunicator::getMPIType( stype );

    SCAI_MPICALL( logger,
        MPI_File_read_at_all( mFileHandle, offset, val, n, commType, &stat ),
        "MPI_File_read_at_all" )

    int readEntries = 0;

    SCAI_MPICALL( logger, MPI_Get_count( &stat, commType, &readEntries ),
                  "MPI_Get_count<" << stype << ">" )

    return static_cast<size_t>( readEntries );
}

/* -------------------------------------------------------------------------- */

size_t MPICollectiveFile::getSize() const
{
    MPI_Offset size;
    SCAI_MPICALL( logger, MPI_File_get_size( mFileHandle, &size ), "MPI_File_get_size" )
    return static_cast<size_t>( size );
}

/* -------------------------------------------------------------------------- */

}

}
