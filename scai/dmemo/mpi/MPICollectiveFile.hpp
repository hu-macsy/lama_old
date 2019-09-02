/**
 * @file MPICollectiveFile.hpp
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

#pragma once

#include <mpi.h>

#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/dmemo/mpi/MPICommunicator.hpp>
#include <scai/hmemo.hpp>
#include <scai/hmemo/HostWriteOnlyAccess.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/logging.hpp>

namespace scai
{

namespace dmemo
{

/**
 *   MPI implementation of a collective file using MPI parallel I/O routines.
 */
class COMMON_DLL_IMPORTEXPORT MPICollectiveFile : public CollectiveFile

{
public:

    /**
     *   Constructor of a file
     */
    MPICollectiveFile( CommunicatorPtr comm );

    /**
     *   Implementation of pure method CollectiveFile::open
     */
    virtual void open( const char* fileName, const char* fileMode );

    /**
     *   Implementation of pure method CollectiveFile::close
     */
    virtual void close();

    /**
     *   Implementation of pure method CollectiveFile::writeSingleImpl
     */
    virtual size_t writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype );

    /**
     *   Implementation of pure method CollectiveFile::writeAllImpl
     */
    virtual size_t writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype );

    /**
     *   Implementation of pure method CollectiveFile::readSingleImpl
     */
    virtual size_t readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype );

    /**
     *   Implementation of pure method CollectiveFile::readAllImpl
     */
    virtual size_t readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype );

    /**
     *   Return the size of the file (in number of bytes)
     */
    virtual size_t getSize() const;

private:

    MPI_File mFileHandle;
};

}

}
