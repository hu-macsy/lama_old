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
    // we need the MPI communicator from the communicator

    auto mpiComm = std::dynamic_pointer_cast<const MPICommunicator>( mComm );

    SCAI_ASSERT( mpiComm, "MPICollective file only with MPI communicator, comm = " << *mComm )

    int openMode = MPI_MODE_WRONLY | MPI_MODE_CREATE;

    if ( strcmp( fileMode, "r" ) == 0 )
    {
        openMode = MPI_MODE_RDONLY;
    }

    SCAI_MPICALL( logger,
                  MPI_File_open( mpiComm->getMPIComm(), const_cast<char*>( fileName ), openMode, MPI_INFO_NULL, &mFileHandle ),
                  "MPI_File_open" )

    CollectiveFile::set( fileName, 0 );
}

/* -------------------------------------------------------------------------- */

void MPICollectiveFile::close()
{
    SCAI_LOG_INFO( logger, "close file " << mFileName << ", offset = " << mOffset );

    SCAI_MPICALL( logger, MPI_File_close( &mFileHandle ), "MPI_File_close" )
}

/* -------------------------------------------------------------------------- */

void MPICollectiveFile::writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    MPI_Status stat;
    SCAI_MPICALL( logger, 
      MPI_File_write_at( mFileHandle, offset, const_cast<void*>( val ), n, MPICommunicator::getMPIType( stype ), &stat ),
      "MPI_File_write_at" )
}

void MPICollectiveFile::writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype )
{
    MPI_Status stat;

    SCAI_MPICALL( logger,
        MPI_File_write_at_all(
            mFileHandle, 
            offset, 
            const_cast<void*>( val ), n,
            MPICommunicator::getMPIType( stype ), &stat ),
        "MPI_File_write_at_all" )
}

/* -------------------------------------------------------------------------- */

void MPICollectiveFile::readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    MPI_Status stat;
    SCAI_MPICALL( logger,
        MPI_File_read_at( mFileHandle, offset, val, n, MPICommunicator::getMPIType( stype ), &stat ),
        "MPI_File_read_at" )
}

void MPICollectiveFile::readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype )
{
    MPI_Status stat;
    SCAI_MPICALL( logger,
        MPI_File_read_at_all( mFileHandle, offset, val, n, MPICommunicator::getMPIType( stype ), &stat ),
        "MPI_File_read_at_all" )
}

/* -------------------------------------------------------------------------- */

}

}
