/**
 * @file NoCollectiveFile.hpp
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

#include <scai/dmemo/CollectiveFile.hpp>

#include <scai/logging.hpp>

#include <cstdio>

namespace scai
{

namespace dmemo
{

class COMMON_DLL_IMPORTEXPORT NoCollectiveFile : public CollectiveFile

{
public:

    NoCollectiveFile( CommunicatorPtr comm );

    virtual ~NoCollectiveFile();

    virtual void open( const char* fileName, const char* fileMode );

    virtual void close();

    virtual void writeSingleImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype );

    virtual void writeAllImpl( const size_t offset, const void* val, const size_t n, const common::ScalarType stype );

    virtual void readSingleImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype );

    virtual void readAllImpl( void* val, const size_t n, const size_t offset, const common::ScalarType stype );

private:

    FILE* mFile;
    
};

}

}
